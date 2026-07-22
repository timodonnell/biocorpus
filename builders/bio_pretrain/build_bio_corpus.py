#!/usr/bin/env python3
# Copyright (c) 2026. Released under the Apache License 2.0.
"""
Proof-of-concept builder for a *sequence-first* biology pre-training corpus
from UniProtKB/Swiss-Prot and Ensembl.

Design rationale (see ../../analyses/thebiocollection-review.html for the review
that motivated these choices):

  * SEQUENCE-FIRST ordering. Each document is rendered as
        >{db}:{accession} {name} [{organism}]
        <protein> ...sequence... </protein>

        <clean natural-language metadata about the sequence>
    In an autoregressive LM the loss factorises left-to-right, so putting the
    sequence first trains P(metadata | sequence) -- the *recognition* direction
    ("given a sequence, tell me what it is"), which is the useful inference
    direction for a science assistant. `--ordering metadata_first` trains the
    reverse (design) direction; `--ordering both` emits one of each.

  * PROVENANCE is a first-class column (source, version, url, license,
    accession, organism, taxid). TheBioCollection shipped only {text,
    record_type}; that made filtering, dedup, ablation and licensing
    impossible. We keep the structured fields *and* the rendered `text`.

  * NO EMPTY-SLOT BUGS. Missing fields are omitted, never rendered as blanks.
    (~5% of TheBioCollection records contained "performed at pH  and .")

  * DEDUP + CAPS. Exact-sequence dedup (a stand-in for UniRef clustering) and
    a hard record cap keep the pilot small and non-redundant.

This is a POC: exact dedup only (not UniRef50/90), no distributed execution,
and it streams the *head* of each source. It is meant to make the schema and
the sequence-first rendering concrete and runnable, not to be a full pipeline.

Usage:
    # rich, curated protein records fetched individually via the UniProt API
    python build_bio_corpus.py uniprot --accessions P02768,P00533,P69905,P0DTC2 \
        --ordering sequence_first --out samples/uniprot_swissprot.sample.jsonl

    # stream the head of the full Swiss-Prot flat file instead
    python build_bio_corpus.py uniprot --limit 200 --out out.jsonl

    # Ensembl peptides (default: yeast); use --seq-type cdna for nucleotide transcripts
    python build_bio_corpus.py ensembl --limit 40 --out samples/ensembl_pep.sample.jsonl
"""
from __future__ import annotations

import argparse
import dataclasses
import gzip
import hashlib
import io
import json
import os
import re
import sys
import urllib.request
from dataclasses import dataclass, field
from typing import Iterator, Optional

# --------------------------------------------------------------------------- #
# Schema                                                                       #
# --------------------------------------------------------------------------- #

SCHEMA_VERSION = "0.1"

SEQ_DELIMS = {"aa": ("<protein>", "</protein>"), "dna": ("<dna>", "</dna>")}

LICENSES = {
    "uniprot_swissprot": "CC-BY-4.0 (https://www.uniprot.org/help/license)",
    # Ensembl's own annotation is unrestricted; some genome assemblies carry
    # third-party terms, so downstream users should verify per species.
    "ensembl": "Open, no restrictions on Ensembl annotation "
    "(https://www.ensembl.org/info/about/legal/disclaimer.html); verify assembly terms per species",
}


@dataclass
class BioRecord:
    """One sequence entry + its metadata + the rendered pre-training document.

    The `text` field is what a tokenizer would consume. Every other field is
    provenance/metadata that TheBioCollection dropped and that we deliberately
    keep so the corpus can be filtered, deduped, ablated and license-checked.
    """

    id: str  # globally unique, e.g. "uniprot:P02768"
    source: str  # "uniprot_swissprot" | "ensembl"
    source_version: str  # release label
    source_url: str  # exact URL the entry came from
    license: str
    accession: str
    entity_type: str  # "protein" | "transcript"
    seq_type: str  # "aa" | "dna"
    seq_len: int
    organism: Optional[str]
    taxid: Optional[str]
    gene: Optional[str]
    name: Optional[str]
    annotations: dict  # structured, source-specific
    sequence: str
    ordering: str = "sequence_first"
    schema_version: str = SCHEMA_VERSION
    text: str = ""  # filled by render()


# --------------------------------------------------------------------------- #
# Rendering (shared by both sources)                                           #
# --------------------------------------------------------------------------- #

_EVID = re.compile(r"\s*\{ECO:[^}]*\}")


def _clean(s: Optional[str]) -> Optional[str]:
    """Strip UniProt evidence tags and squeeze whitespace."""
    if not s:
        return None
    s = _EVID.sub("", s).strip().rstrip(".").strip()
    s = re.sub(r"\s+", " ", s)
    return s or None


def _anchor(rec: BioRecord) -> str:
    prefix = "sp" if rec.source == "uniprot_swissprot" else "ensembl"
    bits = [f">{prefix}:{rec.accession}"]
    if rec.name:
        bits.append(rec.name)
    if rec.organism:
        bits.append(f"[{rec.organism}]")
    return " ".join(bits)


def _seq_block(rec: BioRecord) -> str:
    o, c = SEQ_DELIMS[rec.seq_type]
    return f"{o}{rec.sequence}{c}"


def _metadata_lines(rec: BioRecord) -> list[str]:
    """Order-neutral, clean natural-language metadata. Missing fields omitted."""
    a = rec.annotations
    unit = "residue" if rec.seq_type == "aa" else "nt"
    lines: list[str] = []

    # Lead identity sentence.
    ident = rec.name or rec.accession
    who = f"{ident} — {_db_label(rec)} {rec.accession} — is a {rec.seq_len}-{unit} "
    who += "protein" if rec.entity_type == "protein" else "transcript"
    if rec.organism:
        who += f" from {rec.organism}"
        if rec.taxid:
            who += f" (NCBI taxon {rec.taxid})"
    if rec.gene:
        who += f"; gene {rec.gene}"
    lines.append(who + ".")

    for label, key in (
        ("Function", "function"),
        ("Catalytic activity", "catalytic_activity"),
        ("Subcellular location", "subcellular_location"),
        ("Pathway", "pathway"),
        ("Involvement in disease", "disease"),
        ("Description", "description"),
        ("Location", "location"),
        ("Biotype", "biotype"),
    ):
        v = a.get(key)
        if v:
            lines.append(f"{label}: {v}")

    if a.get("go"):
        lines.append("GO annotations: " + a["go"])
    if a.get("features"):
        lines.append("Sequence features: " + a["features"])
    if a.get("keywords"):
        lines.append("Keywords: " + a["keywords"])
    if a.get("xrefs"):
        lines.append("Cross-references: " + a["xrefs"])
    if a.get("lineage"):
        lines.append("Lineage: " + a["lineage"])
    return lines


def _db_label(rec: BioRecord) -> str:
    return {
        "uniprot_swissprot": "UniProtKB/Swiss-Prot",
        "ensembl": "Ensembl",
    }.get(rec.source, rec.source)


def render(rec: BioRecord, ordering: str) -> BioRecord:
    """Populate rec.text according to `ordering`; returns rec."""
    anchor, seq, meta = _anchor(rec), _seq_block(rec), "\n".join(_metadata_lines(rec))
    if ordering == "metadata_first":
        rec.text = f"{meta}\n\n{anchor}\n{seq}"
    else:  # sequence_first (default)
        rec.text = f"{anchor}\n{seq}\n\n{meta}"
    rec.ordering = ordering
    return rec


# --------------------------------------------------------------------------- #
# IO helpers                                                                   #
# --------------------------------------------------------------------------- #

_UA = {"User-Agent": "biocorpus-poc/0.1 (research)"}


def open_text(path_or_url: str) -> io.TextIOBase:
    """Return a streaming text handle for a local or remote (optionally gz) source."""
    if path_or_url.startswith(("http://", "https://", "ftp://")):
        raw = urllib.request.urlopen(urllib.request.Request(path_or_url, headers=_UA), timeout=120)
        stream: io.BufferedReader = raw  # type: ignore[assignment]
    else:
        stream = open(path_or_url, "rb")
    if path_or_url.endswith(".gz"):
        stream = gzip.GzipFile(fileobj=stream)  # type: ignore[assignment]
    return io.TextIOWrapper(stream, encoding="ascii", errors="replace")


def fetch_text(url: str) -> str:
    return urllib.request.urlopen(urllib.request.Request(url, headers=_UA), timeout=120).read().decode(
        "ascii", "replace"
    )


# --------------------------------------------------------------------------- #
# Source: UniProtKB/Swiss-Prot                                                 #
# --------------------------------------------------------------------------- #

UNIPROT_DAT_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot_sprot.dat.gz"
)
_GO_ASPECT = {"C": "component", "F": "function", "P": "process"}
_FEATURE_INTEREST = (
    "SIGNAL", "TRANSMEM", "DOMAIN", "ACT_SITE", "BINDING", "METAL",
    "SITE", "MOD_RES", "DISULFID", "CARBOHYD", "MOTIF", "REGION",
)


def _gene_from_uniprot(gene_name) -> Optional[str]:
    if not gene_name:
        return None
    g = gene_name[0] if isinstance(gene_name, list) else gene_name
    if isinstance(g, dict):
        for k in ("Name", "OrderedLocusNames", "ORFNames", "Synonyms"):
            v = g.get(k)
            if v:
                return _clean(v[0] if isinstance(v, list) else v)
    return _clean(g) if isinstance(g, str) else None


def _summarise_features(features) -> Optional[str]:
    from collections import Counter, defaultdict

    counts: Counter = Counter()
    examples: dict[str, list[str]] = defaultdict(list)
    for ft in features or []:
        t = getattr(ft, "type", None)
        if not t:
            continue
        counts[t] += 1
        if t in _FEATURE_INTEREST and len(examples[t]) < 2:
            try:
                start = int(ft.location.start) + 1
                end = int(ft.location.end)
                span = f"{start}" if start == end else f"{start}-{end}"
            except Exception:
                span = "?"
            note = None
            q = getattr(ft, "qualifiers", None)
            if isinstance(q, dict):
                note = _clean(q.get("note"))
            examples[t].append(f"{span}" + (f" ({note})" if note else ""))
    if not counts:
        return None
    parts = []
    for t, n in counts.most_common(8):
        label = t.lower().replace("_", " ")
        piece = f"{n} {label}"
        if examples.get(t):
            piece += " [" + "; ".join(examples[t]) + "]"
        parts.append(piece)
    return "; ".join(parts)


def _uniprot_record_from_swiss(rec, source_url: str, version: str) -> BioRecord:
    acc = rec.accessions[0]
    name = None
    m = re.search(r"Full=([^;{]+)", rec.description or "")
    if m:
        name = _clean(m.group(1))

    comments: dict[str, str] = {}
    for c in rec.comments or []:
        if ": " in c:
            topic, body = c.split(": ", 1)
            comments.setdefault(topic.strip(), _clean(body) or "")

    go_terms = []
    xref_by_db: dict[str, list[str]] = {}
    for x in rec.cross_references or []:
        db = x[0]
        if db == "GO" and len(x) >= 3:
            aspect, _, term = x[2].partition(":")
            go_terms.append(f"{term} ({_GO_ASPECT.get(aspect, aspect)})")
        elif db in ("PDB", "Pfam", "InterPro", "KEGG", "Reactome", "RefSeq", "EMBL"):
            xref_by_db.setdefault(db, []).append(x[1])
    xref_txt = "; ".join(
        f"{db} " + ", ".join(ids[:4]) + (f" +{len(ids) - 4} more" if len(ids) > 4 else "")
        for db, ids in xref_by_db.items()
    )

    annotations = {
        "function": comments.get("FUNCTION"),
        "catalytic_activity": comments.get("CATALYTIC ACTIVITY"),
        "subcellular_location": comments.get("SUBCELLULAR LOCATION"),
        "pathway": comments.get("PATHWAY"),
        # UniProt DISEASE comments append a "Note=..." clause; keep only the definition.
        "disease": _clean(comments["DISEASE"].split(" Note=")[0]) if comments.get("DISEASE") else None,
        "go": "; ".join(go_terms[:12]) + (f"; +{len(go_terms) - 12} more" if len(go_terms) > 12 else "")
        if go_terms
        else None,
        "features": _summarise_features(rec.features),
        "keywords": "; ".join(rec.keywords[:12]) if rec.keywords else None,
        "xrefs": xref_txt or None,
        "lineage": " > ".join((rec.organism_classification or [])[:8]) or None,
        "entry_name": rec.entry_name,
    }
    annotations = {k: v for k, v in annotations.items() if v}

    return BioRecord(
        id=f"uniprot:{acc}",
        source="uniprot_swissprot",
        source_version=version,
        source_url=source_url,
        license=LICENSES["uniprot_swissprot"],
        accession=acc,
        entity_type="protein",
        seq_type="aa",
        seq_len=rec.sequence_length,
        organism=_clean(rec.organism),
        taxid=(rec.taxonomy_id or [None])[0],
        gene=_gene_from_uniprot(getattr(rec, "gene_name", None)),
        name=name,
        annotations=annotations,
        sequence=rec.sequence,
    )


def iter_uniprot(args) -> Iterator[BioRecord]:
    from Bio import SwissProt

    if args.accessions:
        for acc in [a.strip() for a in args.accessions.split(",") if a.strip()]:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.txt"
            handle = io.StringIO(fetch_text(url))
            for rec in SwissProt.parse(handle):
                yield _uniprot_record_from_swiss(rec, url, args.version or "uniprot_rest")
    else:
        src = args.dat or UNIPROT_DAT_URL
        handle = open_text(src)
        for rec in SwissProt.parse(handle):
            yield _uniprot_record_from_swiss(rec, src, args.version or "current_release")


# --------------------------------------------------------------------------- #
# Source: Ensembl (FASTA with rich headers)                                    #
# --------------------------------------------------------------------------- #

ENSEMBL_DEFAULT_PEP = (
    "https://ftp.ensembl.org/pub/release-112/fasta/saccharomyces_cerevisiae/"
    "pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
)
_ENS_KV = re.compile(r"(\w+):(\S+)")


def _parse_ensembl_header(header: str) -> dict:
    """Parse an Ensembl FASTA header into fields.

    Example:
      YLL050C pep chromosome:R64-1-1:XII:39804:40414:-1 gene:YLL050C
      transcript:YLL050C_mRNA gene_biotype:protein_coding
      transcript_biotype:protein_coding gene_symbol:COF1 description:Cofilin, ... [Source:SGD;Acc:...]
    """
    header = header[1:].strip() if header.startswith(">") else header.strip()
    out: dict[str, str] = {}
    desc = ""
    if " description:" in header:
        header, desc = header.split(" description:", 1)
    toks = header.split()
    out["id"] = toks[0]
    out["molecule"] = toks[1] if len(toks) > 1 else ""
    for tok in toks[2:]:
        if ":" in tok:
            k, _, v = tok.partition(":")
            out[k] = v
    if "chromosome" in out:  # assembly:seqname:start:end:strand
        parts = out["chromosome"].split(":")
        if len(parts) >= 5:
            out["assembly"], out["seqname"] = parts[0], parts[1]
            out["location"] = f"{parts[1]}:{parts[2]}-{parts[3]} ({'+' if parts[4] == '1' else '-'})"
    if desc:
        out["description"] = desc.strip()
    return out


def _version_from_url(url: str) -> str:
    m = re.search(r"release-(\d+)", url)
    return f"release-{m.group(1)}" if m else "ensembl"


def iter_ensembl(args) -> Iterator[BioRecord]:
    url = args.url or ENSEMBL_DEFAULT_PEP
    seq_type = "aa" if args.seq_type == "pep" else "dna"
    entity = "protein" if args.seq_type == "pep" else "transcript"
    version = _version_from_url(url)
    handle = open_text(url)

    hdr: Optional[dict] = None
    seq: list[str] = []

    def _emit(h, s):
        sequence = "".join(s)
        acc = h.get("id", "?")
        desc = _clean(h.get("description"))
        if desc:  # the [Source:...] tag is provenance noise in prose; drop it from the description
            desc = re.sub(r"\s*\[Source:[^\]]*\]", "", desc).strip()
        annotations = {
            "description": desc,
            "biotype": h.get("transcript_biotype") or h.get("gene_biotype"),
            "location": h.get("location"),
            "gene_id": h.get("gene"),
            "transcript_id": h.get("transcript"),
            "assembly": h.get("assembly"),
        }
        annotations = {k: v for k, v in annotations.items() if v}
        return BioRecord(
            id=f"ensembl:{acc}",
            source="ensembl",
            source_version=version,
            source_url=url,
            license=LICENSES["ensembl"],
            accession=acc,
            entity_type=entity,
            seq_type=seq_type,
            seq_len=len(sequence),
            organism=args.organism,
            taxid=args.taxid,
            gene=h.get("gene_symbol") or h.get("gene"),
            name=h.get("gene_symbol"),  # short identity for the anchor; full text is in annotations.description
            annotations=annotations,
            sequence=sequence,
        )

    for line in handle:
        if line.startswith(">"):
            if hdr is not None:
                yield _emit(hdr, seq)
            hdr, seq = _parse_ensembl_header(line), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        yield _emit(hdr, seq)


# --------------------------------------------------------------------------- #
# Build driver                                                                 #
# --------------------------------------------------------------------------- #


def build(args) -> None:
    it = iter_uniprot(args) if args.source == "uniprot" else iter_ensembl(args)

    seen_seq: set[str] = set()
    n_written = n_dup = n_toolong = n_in = 0
    orderings = ["sequence_first", "metadata_first"] if args.ordering == "both" else [args.ordering]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        for rec in it:
            n_in += 1
            if args.max_seq_len and rec.seq_len > args.max_seq_len:
                n_toolong += 1
                continue
            h = hashlib.sha1(rec.sequence.encode()).hexdigest()
            if h in seen_seq:
                n_dup += 1
                continue
            seen_seq.add(h)
            for i, ordering in enumerate(orderings):
                r = render(dataclasses.replace(rec), ordering)
                if len(orderings) > 1:
                    r.id = f"{rec.id}#{ordering[:2]}"  # keep ids unique
                fh.write(json.dumps(dataclasses.asdict(r), ensure_ascii=False) + "\n")
            n_written += 1
            if args.limit and n_written >= args.limit:
                break

    print(
        f"[{args.source}] wrote {n_written} entries "
        f"({n_written * len(orderings)} docs, ordering={args.ordering}) -> {args.out}\n"
        f"          read={n_in} dup_seq_skipped={n_dup} too_long_skipped={n_toolong}",
        file=sys.stderr,
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="source", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", required=True, help="output .jsonl path")
    common.add_argument("--limit", type=int, default=50, help="max unique entries to write")
    common.add_argument(
        "--ordering",
        choices=["sequence_first", "metadata_first", "both"],
        default="sequence_first",
    )
    common.add_argument("--max-seq-len", type=int, default=0, help="skip sequences longer than this (0=off)")

    up = sub.add_parser("uniprot", parents=[common], help="UniProtKB/Swiss-Prot")
    up.add_argument("--accessions", help="comma-separated accessions (fetched via UniProt REST)")
    up.add_argument("--dat", help="local/remote uniprot_sprot.dat[.gz] (default: UniProt FTP head)")
    up.add_argument("--version", help="source_version label")

    en = sub.add_parser("ensembl", parents=[common], help="Ensembl FASTA (rich headers)")
    en.add_argument("--url", help="Ensembl pep/cdna FASTA .gz URL (default: yeast pep)")
    en.add_argument("--seq-type", choices=["pep", "cdna"], default="pep")
    en.add_argument("--organism", default="Saccharomyces cerevisiae")
    en.add_argument("--taxid", default="4932")

    args = p.parse_args()
    build(args)


if __name__ == "__main__":
    main()
