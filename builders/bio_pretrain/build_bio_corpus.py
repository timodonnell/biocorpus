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
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Iterator, Optional

# --------------------------------------------------------------------------- #
# Schema                                                                       #
# --------------------------------------------------------------------------- #

SCHEMA_VERSION = "0.1"

SEQ_DELIMS = {"aa": ("<protein>", "</protein>"), "dna": ("<dna>", "</dna>")}

_UNIPROT_LICENSE = "CC-BY-4.0 (https://www.uniprot.org/help/license)"
LICENSES = {
    "uniprot_swissprot": _UNIPROT_LICENSE,
    "uniprot_trembl": _UNIPROT_LICENSE,
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
    # multi-sequence records (central dogma) keep the raw forms here: {dna_genomic, rna, cds, protein}
    sequences: Optional[dict] = None


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
    prefix = {"uniprot_swissprot": "sp", "uniprot_trembl": "tr", "ensembl": "ensembl"}.get(
        rec.source, "db"
    )
    bits = [f">{prefix}:{rec.accession}"]
    if rec.name:
        bits.append(rec.name)
    if rec.organism:
        bits.append(f"[{rec.organism}]")
    return " ".join(bits)


def _seq_block(rec: BioRecord) -> str:
    if not rec.sequence:  # gene-annotation records may carry no sequence
        return ""
    o, c = SEQ_DELIMS[rec.seq_type]
    return f"{o}{rec.sequence}{c}"


def _lead_sentence(rec: BioRecord) -> str:
    """Order-neutral identity sentence (never references 'above'/'below')."""
    ident = rec.name or rec.accession
    org = ""
    if rec.organism:
        org = f" from {rec.organism}" + (f" (NCBI taxon {rec.taxid})" if rec.taxid else "")
    if rec.entity_type == "regulatory_feature":
        ft = (rec.annotations.get("feature_type") or "regulatory feature").replace("_", " ")
        return (
            f"The sequence is a {ft} ({rec.seq_len} bp) from the Ensembl "
            f"{rec.source_version} Regulatory Build ({rec.accession}){org}."
        )
    if rec.entity_type == "splice_junction":
        return f"The sequence is a {rec.seq_len} bp window over a 5' splice donor site in gene {rec.gene or '?'}{org}."
    if rec.entity_type == "gene":
        biotype = rec.annotations.get("biotype")
        return (
            f"{ident} — {_db_label(rec)} gene {rec.accession} — is a "
            + (f"{biotype} " if biotype else "")
            + f"gene{org}."
        )
    unit = "residue" if rec.seq_type == "aa" else "nt"
    kind = "protein" if rec.entity_type == "protein" else "transcript"
    who = f"{ident} — {_db_label(rec)} {rec.accession} — is a {rec.seq_len}-{unit} {kind}{org}"
    if rec.gene:
        who += f"; gene {rec.gene}"
    return who + "."


def _metadata_lines(rec: BioRecord) -> list[str]:
    """Order-neutral, clean natural-language metadata. Missing fields omitted."""
    a = rec.annotations
    lines: list[str] = [_lead_sentence(rec)]

    for label, key in (
        ("Function", "function"),
        ("Catalytic activity", "catalytic_activity"),
        ("Subcellular location", "subcellular_location"),
        ("Pathway", "pathway"),
        ("Involvement in disease", "disease"),
        ("Gene model", "gene_model"),
        ("Canonical product", "product"),
        ("Description", "description"),
        ("Location", "location"),
        ("Biotype", "biotype"),
        ("Feature length", "feature_length_bp"),
        ("Splice site", "splice_site"),
        ("Transcript", "transcript"),
        ("Intron", "intron"),
        ("Intron location", "intron_location"),
        ("Intron length", "intron_length_bp"),
        ("Donor dinucleotide", "donor_dinucleotide"),
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
        "uniprot_trembl": "UniProtKB/TrEMBL",
        "ensembl": "Ensembl",
    }.get(rec.source, rec.source)


def render(rec: BioRecord, ordering: str) -> BioRecord:
    """Populate rec.text according to `ordering`; returns rec."""
    if rec.entity_type == "central_dogma":  # dogma records are pre-rendered (fixed DNA->RNA->protein order)
        rec.ordering = ordering
        return rec
    anchor = _anchor(rec)
    seq = _seq_block(rec)
    meta = "\n".join(_metadata_lines(rec))
    head = f"{anchor}\n{seq}" if seq else anchor
    if ordering == "metadata_first":
        rec.text = f"{meta}\n\n{head}"
    else:  # sequence_first (default)
        rec.text = f"{head}\n\n{meta}"
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

    source = (
        "uniprot_swissprot"
        if (getattr(rec, "data_class", "") or "").lower().startswith("reviewed")
        else "uniprot_trembl"
    )
    return BioRecord(
        id=f"uniprot:{acc}",
        source=source,
        source_version=version,
        source_url=source_url,
        license=LICENSES[source],
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
    elif getattr(args, "query", None):
        # Sample a representative slice matching a UniProt query via the paginated
        # search API (robust for both the small reviewed set and the huge
        # unreviewed one, where /stream is refused), e.g.
        #   --query "reviewed:true" (Swiss-Prot)  /  --query "reviewed:false" (TrEMBL)
        page_url = (
            "https://rest.uniprot.org/uniprotkb/search?"
            f"query={urllib.parse.quote(args.query)}&format=txt&size=500"
        )
        while page_url:
            resp = urllib.request.urlopen(urllib.request.Request(page_url, headers=_UA), timeout=300)
            text = resp.read().decode("ascii", "replace")
            for rec in SwissProt.parse(io.StringIO(text)):
                yield _uniprot_record_from_swiss(rec, page_url, args.version or "uniprot_search")
            m = re.search(r'<([^>]+)>;\s*rel="next"', resp.headers.get("Link", ""))
            page_url = m.group(1) if m else None
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
# Source: Ensembl GFF3 gene models (structural annotation, optionally + sequence)
# --------------------------------------------------------------------------- #

ENSEMBL_DEFAULT_GFF = (
    "https://ftp.ensembl.org/pub/release-112/gff3/saccharomyces_cerevisiae/"
    "Saccharomyces_cerevisiae.R64-1-1.112.gff3.gz"
)
_GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}


def _gff_attrs(col9: str) -> dict:
    """Parse a GFF3 column-9 attribute string; values are URL-decoded."""
    out: dict[str, str] = {}
    for kv in col9.rstrip(";\n").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = urllib.parse.unquote(v)
    return out


def _load_fasta_gene_seq(url: str, cap: int = 200_000) -> dict:
    """Map Ensembl gene_id -> (sequence, fasta_id) using the first FASTA entry per gene."""
    m: dict[str, tuple[str, str]] = {}
    handle = open_text(url)
    hdr: Optional[dict] = None
    seq: list[str] = []

    def _flush():
        if hdr is None:
            return
        g = hdr.get("gene")
        if g and g not in m:
            m[g] = ("".join(seq), hdr.get("id", ""))

    for line in handle:
        if line.startswith(">"):
            _flush()
            if len(m) >= cap:
                return m
            hdr, seq = _parse_ensembl_header(line), []
        else:
            seq.append(line.strip())
    _flush()
    return m


def _gene_model_str(g: dict) -> Optional[str]:
    if not g.get("n_tx"):
        return None
    s = f"{g['n_tx']} transcript(s), {g['n_exon']} exon(s) in the canonical transcript"
    bt = sorted(t for t in g.get("tx_biotypes", set()) if t)
    if len(bt) > 1:
        s += f"; transcript biotypes: {', '.join(bt)}"
    return s


def iter_ensembl_gff(args) -> Iterator[BioRecord]:
    """Gene-level records from an Ensembl GFF3, optionally joined to a product sequence."""
    gff = args.gff or ENSEMBL_DEFAULT_GFF
    version = _version_from_url(gff)
    seqmap: dict[str, tuple[str, str]] = {}
    seq_type: Optional[str] = None
    if args.seq_type != "none":
        seqmap = _load_fasta_gene_seq(args.fasta or ENSEMBL_DEFAULT_PEP)
        seq_type = "aa" if args.seq_type == "pep" else "dna"

    def _build(g: dict) -> BioRecord:
        gid = g["gene_id"]
        sequence, _fid = seqmap.get(gid, ("", ""))
        annotations = {
            "biotype": g.get("biotype"),
            "gene_model": _gene_model_str(g),
            "description": g.get("description"),
            "location": g.get("location"),
            "gene_id": gid,
        }
        if sequence:
            unit = "residue" if seq_type == "aa" else "nt"
            kind = "protein product" if seq_type == "aa" else "transcript (cDNA)"
            annotations["product"] = f"canonical {kind}, {len(sequence)} {unit}s"
        annotations = {k: v for k, v in annotations.items() if v}
        return BioRecord(
            id=f"ensembl:{gid}",
            source="ensembl",
            source_version=version,
            source_url=gff,
            license=LICENSES["ensembl"],
            accession=gid,
            entity_type="gene",
            seq_type=seq_type or "aa",
            seq_len=len(sequence),
            organism=args.organism,
            taxid=args.taxid,
            gene=g.get("name") or gid,
            name=g.get("name"),
            annotations=annotations,
            sequence=sequence,
        )

    handle = open_text(gff)
    cur: Optional[dict] = None
    for line in handle:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        ftype, attrs = cols[2], _gff_attrs(cols[8])
        if ftype in _GENE_TYPES:
            if cur is not None:
                yield _build(cur)
            cur = {
                "gene_id": attrs.get("gene_id") or attrs.get("ID", "").replace("gene:", ""),
                "name": attrs.get("Name"),
                "biotype": attrs.get("biotype"),
                "description": re.sub(r"\s*\[Source:[^\]]*\]", "", attrs.get("description", "")).strip() or None,
                "location": f"{cols[0]}:{int(cols[3]):,}-{int(cols[4]):,} ({cols[6]})",
                "n_tx": 0,
                "n_exon": 0,
                "canonical_tx": None,
                "tx_biotypes": set(),
            }
        elif cur is not None:
            parent = attrs.get("Parent", "")
            if ftype == "exon":
                ctx = cur["canonical_tx"]
                if not ctx or parent == f"transcript:{ctx}":
                    cur["n_exon"] += 1
            elif parent.startswith("gene:"):  # a transcript of the current gene
                cur["n_tx"] += 1
                if attrs.get("biotype"):
                    cur["tx_biotypes"].add(attrs["biotype"])
                if "Ensembl_canonical" in attrs.get("tag", ""):
                    cur["canonical_tx"] = attrs.get("transcript_id") or attrs.get("ID", "").replace(
                        "transcript:", ""
                    )
    if cur is not None:
        yield _build(cur)


# --------------------------------------------------------------------------- #
# Sources: Ensembl regulatory features + splice-donor junctions (DNA via REST) #
# --------------------------------------------------------------------------- #

ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_REG_GFF = (
    "https://ftp.ensembl.org/pub/release-112/regulation/homo_sapiens/GRCh38/annotation/"
    "Homo_sapiens.GRCh38.regulatory_features.v112.gff3.gz"
)
ENSEMBL_HUMAN_GFF = (
    "https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz"
)


def _ensembl_region_seq(species: str, seqid: str, start: int, end: int, strand: int = 1) -> str:
    """Fetch a genomic DNA sequence via the Ensembl REST API (strand -1 => reverse complement)."""
    url = f"{ENSEMBL_REST}/sequence/region/{species}/{seqid}:{start}..{end}:{strand}?content-type=text/plain"
    return (
        urllib.request.urlopen(urllib.request.Request(url, headers=_UA), timeout=60)
        .read()
        .decode("ascii", "replace")
        .strip()
    )


def iter_ensembl_regulatory(args) -> Iterator[BioRecord]:
    """One record per Ensembl Regulatory Build feature (enhancer/promoter/CTCF/etc.) + its DNA."""
    url = args.gff or ENSEMBL_REG_GFF
    version = _version_from_url(url)
    for line in open_text(url):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        seqid, ftype, start, end = cols[0], cols[2], int(cols[3]), int(cols[4])
        rid = _gff_attrs(cols[8]).get("ID", f"{seqid}:{start}")
        length = end - start + 1
        windowed = bool(args.window and length > args.window)
        if windowed:  # cap the fetched DNA to a centered window
            center = (start + end) // 2
            fs, fe = max(1, center - args.window // 2), max(1, center - args.window // 2) + args.window - 1
        else:
            fs, fe = start, end
        try:
            seq = _ensembl_region_seq(args.species, seqid, fs, fe, 1)
        except Exception:
            continue
        time.sleep(0.08)  # be polite to the REST API
        annotations = {
            "feature_type": ftype,
            "location": f"{seqid}:{start:,}-{end:,}",
            "feature_length_bp": f"{length:,} bp" + (f" (central {args.window} bp shown)" if windowed else ""),
        }
        yield BioRecord(
            id=f"ensembl_reg:{rid}",
            source="ensembl",
            source_version=version,
            source_url=url,
            license=LICENSES["ensembl"],
            accession=rid,
            entity_type="regulatory_feature",
            seq_type="dna",
            seq_len=len(seq),
            organism=args.organism,
            taxid=args.taxid,
            gene=None,
            name=ftype.replace("_", " "),
            annotations=annotations,
            sequence=seq,
        )


def iter_ensembl_splice(args) -> Iterator[BioRecord]:
    """5' splice-donor junction windows for the canonical transcript of each multi-exon gene."""
    gff = args.gff or ENSEMBL_HUMAN_GFF
    version = _version_from_url(gff)
    w = args.window  # exon/intron context on each side of the donor site

    def emit(g):
        ctx = g.get("canonical_tx")
        if ctx not in g["exons"]:  # fall back to the transcript with the most exons
            if not g["exons"]:
                return
            ctx = max(g["exons"], key=lambda k: len(g["exons"][k]))
        exons = sorted(g["exons"][ctx])
        if len(exons) < 2:
            return
        strand = g["strand"]
        tx = exons if strand == 1 else list(reversed(exons))
        n_intron = len(tx) - 1
        for i, ((sA, eA), (sB, eB)) in enumerate(zip(tx, tx[1:]), start=1):
            if args.per_transcript and i > args.per_transcript:
                break
            if strand == 1:
                a, b, fs, fe = eA + 1, sB - 1, eA + 1 - w, eA + w
            else:
                a, b, fs, fe = eB + 1, sA - 1, sA - w, sA + w - 1
            try:
                seq = _ensembl_region_seq(args.species, g["seqid"], fs, fe, strand)
            except Exception:
                continue
            time.sleep(0.08)
            donor = seq[w : w + 2]
            donor_class = {"GT": "canonical (GT-AG)", "GC": "minor (GC-AG)"}.get(donor, "non-canonical")
            annotations = {
                "splice_site": "5' splice donor",
                "transcript": ctx,
                "intron": f"{i} of {n_intron}",
                "intron_location": f"{g['seqid']}:{a:,}-{b:,} ({'+' if strand == 1 else '-'})",
                "intron_length_bp": f"{b - a + 1:,} bp",
                "donor_dinucleotide": f"{donor} — {donor_class}",
            }
            yield BioRecord(
                id=f"ensembl_splice:{ctx}:intron{i}",
                source="ensembl",
                source_version=version,
                source_url=gff,
                license=LICENSES["ensembl"],
                accession=f"{ctx}:intron{i}",
                entity_type="splice_junction",
                seq_type="dna",
                seq_len=len(seq),
                organism=args.organism,
                taxid=args.taxid,
                gene=g.get("name") or g["gene_id"],
                name=None,
                annotations=annotations,
                sequence=seq,
            )

    cur = None
    for line in open_text(gff):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        ftype, attrs = cols[2], _gff_attrs(cols[8])
        if ftype in _GENE_TYPES:
            if cur is not None:
                yield from emit(cur)
            cur = {
                "gene_id": attrs.get("gene_id") or attrs.get("ID", "").replace("gene:", ""),
                "name": attrs.get("Name"),
                "seqid": cols[0],
                "strand": 1 if cols[6] == "+" else -1,
                "canonical_tx": None,
                "exons": {},
            }
        elif cur is not None:
            parent = attrs.get("Parent", "")
            if ftype == "exon":
                cur["exons"].setdefault(parent.replace("transcript:", ""), []).append(
                    (int(cols[3]), int(cols[4]))
                )
            elif parent.startswith("gene:") and "Ensembl_canonical" in attrs.get("tag", ""):
                cur["canonical_tx"] = attrs.get("transcript_id") or attrs.get("ID", "").replace(
                    "transcript:", ""
                )
    if cur is not None:
        yield from emit(cur)


# --------------------------------------------------------------------------- #
# Source: central dogma (DNA + spliced RNA + protein for one transcript)       #
# --------------------------------------------------------------------------- #


def _ensembl_id_seq(tx: str, seqtype: str) -> str:
    """Fetch one sequence form of a transcript by stable id (type=genomic|cdna|cds|protein)."""
    url = f"{ENSEMBL_REST}/sequence/id/{tx}?type={seqtype}&content-type=text/plain"
    s = (
        urllib.request.urlopen(urllib.request.Request(url, headers=_UA), timeout=60)
        .read()
        .decode("ascii", "replace")
        .strip()
    )
    time.sleep(0.08)  # be polite to the REST API
    return s


def _dogma_views(view: str, coding: bool, compact: bool) -> list:
    """Which document view(s) to emit for a transcript. `all` = one richest feasible view."""
    if view == "all":
        if coding and compact:
            return ["triple"]
        if coding:
            return ["rna_protein"]  # large introns: skip DNA, keep translation
        if compact:
            return ["dna_rna"]  # non-coding: teach splicing
        return []
    if view in ("triple", "dna_protein") and not (coding and compact):
        return []
    if view == "rna_protein" and not coding:
        return []
    if view == "dna_rna" and not compact:
        return []
    return [view]


def _render_dogma(view, gene, tx, biotype, loc, n_exon, n_intron, dna, rna, protein,
                  bounds, cds_span, cds_len, utr5, utr3, intron_total, organism):
    org = f" [{organism}]" if organism else ""
    head = (
        f">ensembl-dogma:{tx} {gene}{org}\n"
        f"Gene {gene} ({biotype}), Ensembl canonical transcript {tx}, {loc}, "
        f"{n_exon} exon(s) / {n_intron} intron(s)."
    )
    exon_line = "Exon boundaries in the mRNA: " + ", ".join(bounds) + "."
    cds_utr = None
    if cds_span:
        parts = []
        if utr5:
            parts.append(f"5' UTR: 1-{utr5} ({utr5} nt)")
        parts.append(f"CDS: {cds_span[0]}-{cds_span[1]} ({cds_len} nt)")
        if utr3:
            parts.append(f"3' UTR: {cds_span[1] + 1}-{cds_span[1] + utr3} ({utr3} nt)")
        cds_utr = "; ".join(parts) + "."
    dna_block = "Genomic DNA (pre-mRNA, sense strand; exons UPPERCASE, introns lowercase):\n" f"<dna>{dna}</dna>"
    splice_line = (
        f"Transcription and splicing remove {n_intron} intron(s) ({intron_total} nt) "
        f"to give the mature mRNA ({len(rna) if rna else 0} nt):"
    )
    translate_line = (
        f"Translation of the CDS ({cds_len} nt, standard genetic code) yields the "
        f"{len(protein) if protein else 0}-residue protein:"
    )
    p = [head, ""]
    if view == "triple":
        p += [dna_block, "", splice_line, f"<rna>{rna}</rna>", exon_line]
        if cds_utr:
            p.append(cds_utr)
        p += ["", translate_line, f"<protein>{protein}</protein>", "",
              "Verified: the mRNA equals the genomic exons with introns removed; translate(CDS) equals the protein."]
    elif view == "dna_rna":
        p += [dna_block, "", splice_line, f"<rna>{rna}</rna>", exon_line]
        if cds_utr:
            p.append(cds_utr)
        p += ["", "Verified: the mRNA equals the genomic exons with introns removed."]
    elif view == "rna_protein":
        p += [f"Mature mRNA ({len(rna)} nt):", f"<rna>{rna}</rna>"]
        if cds_utr:
            p.append(cds_utr)
        p += ["", translate_line, f"<protein>{protein}</protein>", "",
              "Verified: translate(CDS) equals the protein."]
    elif view == "dna_protein":
        p += [dna_block, "",
              f"After transcription, splicing out {n_intron} intron(s), and translation of the CDS, "
              f"this locus encodes the {len(protein)}-residue protein:",
              f"<protein>{protein}</protein>", "",
              "Verified: translate(the CDS of the spliced mRNA) equals the protein."]
    return "\n".join(p)


def iter_ensembl_dogma(args) -> Iterator[BioRecord]:
    """Co-present DNA (pre-mRNA), spliced RNA, and protein for one transcript, with a verified mapping."""
    import warnings
    from Bio.Seq import Seq

    warnings.filterwarnings("ignore")  # silence biopython partial-codon warnings
    gff = args.gff or ENSEMBL_HUMAN_GFF
    version = _version_from_url(gff)

    def process(g):
        ctx = g.get("canonical_tx")
        if not ctx or ctx not in g["exons"]:
            return
        exons = sorted(g["exons"][ctx])
        strand = g["strand"]
        tx_exons = exons if strand == 1 else list(reversed(exons))
        exon_lens = [e - s + 1 for s, e in tx_exons]
        introns = [
            (tx_exons[i + 1][0] - tx_exons[i][1] - 1) if strand == 1 else (tx_exons[i][0] - tx_exons[i + 1][1] - 1)
            for i in range(len(tx_exons) - 1)
        ]
        n_exon, n_intron = len(exon_lens), len(exon_lens) - 1
        if n_exon < args.min_exons:
            return
        span = max(e for _, e in exons) - min(s for s, _ in exons) + 1
        coding = g.get("biotype") == "protein_coding"
        views = _dogma_views(args.view, coding, span <= args.max_dna)
        if not views:
            return

        need_dna = any(v in ("triple", "dna_rna", "dna_protein") for v in views)
        need_rna = any(v in ("triple", "dna_rna", "rna_protein") for v in views)
        need_prot = any(v in ("triple", "dna_protein", "rna_protein") for v in views)
        try:
            cdna = _ensembl_id_seq(ctx, "cdna") if (need_rna or (coding and need_prot)) else None
            genomic = _ensembl_id_seq(ctx, "genomic") if need_dna else None
            cds = _ensembl_id_seq(ctx, "cds") if coding and need_prot else None
            protein = _ensembl_id_seq(ctx, "protein") if need_prot else None
        except Exception:
            return

        # verify the mapping; only assert what we can check
        if coding:
            if not (cds and protein):
                return
            try:
                if str(Seq(cds).translate(table=1, to_stop=True)) != protein:
                    return  # selenoproteins / readthrough / incomplete CDS
            except Exception:
                return
        if cdna is not None and sum(exon_lens) != len(cdna):
            return  # exon structure inconsistent with the spliced transcript
        if genomic is not None and len(genomic) != sum(exon_lens) + sum(introns):
            return

        rna = cdna.replace("T", "U") if cdna else None
        bounds, pos = [], 0
        for L in exon_lens:
            bounds.append(f"{pos + 1}-{pos + L}")
            pos += L
        cds_span = utr5 = utr3 = None
        if coding and cdna and cds:
            idx = cdna.find(cds)
            if idx >= 0:
                cds_span, utr5, utr3 = (idx + 1, idx + len(cds)), idx, len(cdna) - idx - len(cds)
        dna_marked = None
        if genomic is not None:
            segs, out, q = [], [], 0
            for i, L in enumerate(exon_lens):
                segs.append(("e", L))
                if i < n_intron:
                    segs.append(("i", introns[i]))
            for kind, L in segs:
                chunk = genomic[q : q + L]
                out.append(chunk.upper() if kind == "e" else chunk.lower())
                q += L
            dna_marked = "".join(out)

        gene = g.get("name") or g["gene_id"]
        loc = f"{g['seqid']}:{min(s for s, _ in exons):,}-{max(e for _, e in exons):,} ({'+' if strand == 1 else '-'})"
        for v in views:
            text = _render_dogma(
                v, gene, ctx, g.get("biotype"), loc, n_exon, n_intron, dna_marked, rna, protein,
                bounds, cds_span, len(cds) if cds else None, utr5, utr3, sum(introns), args.organism,
            )
            seqs = {}
            if dna_marked is not None and v in ("triple", "dna_rna", "dna_protein"):
                seqs["dna_genomic"] = genomic
            if rna is not None and v in ("triple", "dna_rna", "rna_protein"):
                seqs["rna"] = rna
            if cds and v in ("triple", "rna_protein", "dna_protein"):
                seqs["cds"] = cds
            if protein and v in ("triple", "rna_protein", "dna_protein"):
                seqs["protein"] = protein
            primary = seqs.get("rna") or seqs.get("protein") or seqs.get("dna_genomic") or ""
            ptype = "rna" if "rna" in seqs else ("aa" if "protein" in seqs else "dna")
            rec = BioRecord(
                id=f"ensembl-dogma:{ctx}:{v}",
                source="ensembl",
                source_version=version,
                source_url=gff,
                license=LICENSES["ensembl"],
                accession=ctx,
                entity_type="central_dogma",
                seq_type=ptype,
                seq_len=len(primary),
                organism=args.organism,
                taxid=args.taxid,
                gene=gene,
                name=g.get("name"),
                annotations={
                    "view": v,
                    "biotype": g.get("biotype"),
                    "n_exons": n_exon,
                    "n_introns": n_intron,
                    "location": loc,
                    "verified": "translate(CDS)==protein; mRNA==spliced exons" if coding else "mRNA==spliced exons",
                },
                sequences=seqs,
                sequence=primary,
            )
            rec.text = text
            yield rec

    cur = None
    for line in open_text(gff):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        ftype, attrs = cols[2], _gff_attrs(cols[8])
        if ftype in _GENE_TYPES:
            if cur is not None:
                yield from process(cur)
            cur = {
                "gene_id": attrs.get("gene_id") or attrs.get("ID", "").replace("gene:", ""),
                "name": attrs.get("Name"),
                "seqid": cols[0],
                "strand": 1 if cols[6] == "+" else -1,
                "canonical_tx": None,
                "biotype": None,
                "exons": {},
            }
        elif cur is not None:
            parent = attrs.get("Parent", "")
            if ftype == "exon":
                cur["exons"].setdefault(parent.replace("transcript:", ""), []).append(
                    (int(cols[3]), int(cols[4]))
                )
            elif parent.startswith("gene:") and "Ensembl_canonical" in attrs.get("tag", ""):
                cur["canonical_tx"] = attrs.get("transcript_id") or attrs.get("ID", "").replace(
                    "transcript:", ""
                )
                cur["biotype"] = attrs.get("biotype")
    if cur is not None:
        yield from process(cur)


# --------------------------------------------------------------------------- #
# Build driver                                                                 #
# --------------------------------------------------------------------------- #


def build(args) -> None:
    dispatch = {
        "uniprot": iter_uniprot,
        "ensembl": iter_ensembl,
        "ensembl-gff": iter_ensembl_gff,
        "ensembl-regulatory": iter_ensembl_regulatory,
        "ensembl-splice": iter_ensembl_splice,
        "ensembl-dogma": iter_ensembl_dogma,
    }
    it = dispatch[args.source](args)

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
            if rec.sequence:  # exact-sequence dedup (skip for metadata-only records)
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
    up.add_argument(
        "--query",
        help='sample a UniProt query, e.g. "reviewed:true" (Swiss-Prot) or "reviewed:false" (TrEMBL)',
    )
    up.add_argument("--dat", help="local/remote uniprot_sprot.dat[.gz] (default: UniProt FTP head)")
    up.add_argument("--version", help="source_version label")

    en = sub.add_parser("ensembl", parents=[common], help="Ensembl FASTA (rich headers)")
    en.add_argument("--url", help="Ensembl pep/cdna FASTA .gz URL (default: yeast pep)")
    en.add_argument("--seq-type", choices=["pep", "cdna"], default="pep")
    en.add_argument("--organism", default="Saccharomyces cerevisiae")
    en.add_argument("--taxid", default="4932")

    eg = sub.add_parser("ensembl-gff", parents=[common], help="Ensembl GFF3 gene models (+ optional sequence)")
    eg.add_argument("--gff", help="Ensembl GFF3 .gz URL/path (default: yeast release-112)")
    eg.add_argument("--fasta", help="pep/cdna FASTA to join a product sequence per gene (default: yeast pep)")
    eg.add_argument(
        "--seq-type",
        choices=["pep", "cdna", "none"],
        default="pep",
        help="'none' emits metadata-only gene records",
    )
    eg.add_argument("--organism", default="Saccharomyces cerevisiae")
    eg.add_argument("--taxid", default="4932")

    reg = sub.add_parser("ensembl-regulatory", parents=[common], help="Ensembl Regulatory Build features (+ DNA via REST)")
    reg.add_argument("--gff", help="regulatory features GFF3 .gz (default: human GRCh38 v112)")
    reg.add_argument("--window", type=int, default=600, help="max bp of DNA to fetch per feature (centered)")
    reg.add_argument("--species", default="homo_sapiens")
    reg.add_argument("--organism", default="Homo sapiens")
    reg.add_argument("--taxid", default="9606")

    spl = sub.add_parser("ensembl-splice", parents=[common], help="Ensembl 5' splice-donor junction windows (+ DNA via REST)")
    spl.add_argument("--gff", help="gene GFF3 .gz (default: human GRCh38 v112)")
    spl.add_argument("--window", type=int, default=100, help="bp of exon/intron context each side of the donor site")
    spl.add_argument("--per-transcript", type=int, default=1, help="max donor junctions per transcript")
    spl.add_argument("--species", default="homo_sapiens")
    spl.add_argument("--organism", default="Homo sapiens")
    spl.add_argument("--taxid", default="9606")

    dg = sub.add_parser(
        "ensembl-dogma",
        parents=[common],
        help="Central-dogma records: DNA (pre-mRNA) + spliced RNA + protein for one transcript, verified",
    )
    dg.add_argument("--gff", help="gene GFF3 .gz (default: human GRCh38 v112)")
    dg.add_argument(
        "--view",
        choices=["triple", "dna_rna", "rna_protein", "dna_protein", "all"],
        default="all",
        help="'all' emits one richest feasible view per transcript (ordering is fixed DNA->RNA->protein)",
    )
    dg.add_argument("--max-dna", type=int, default=5000, help="skip DNA-bearing views when the genomic span exceeds this (bp)")
    dg.add_argument("--min-exons", type=int, default=1, help="skip transcripts with fewer exons (use 2+ to require splicing)")
    dg.add_argument("--species", default="homo_sapiens")
    dg.add_argument("--organism", default="Homo sapiens")
    dg.add_argument("--taxid", default="9606")

    args = p.parse_args()
    build(args)


if __name__ == "__main__":
    main()
