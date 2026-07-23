# Sequence-first biology pre-training corpus builder

Turns biological databases into **sequence-first** pre-training documents — raw
sequence, then clean natural-language metadata — with provenance, verified
mappings, and de-duplication. Built to fix the failure modes found in the
[TheBioCollection review](../../analyses/thebiocollection-review.html): no
provenance, ~5% empty-slot template bugs, "verifiable" tasks that aren't
learnable, and DNA/RNA/protein never co-occurring in one document.

One script (`build_bio_corpus.py`), several sources, JSONL out. Only dependency
is `biopython`; stdlib otherwise.

## Sources

| subcommand | what | scale path (billions) |
|---|---|---|
| `uniprot` | Swiss-Prot / TrEMBL proteins (function, GO, features) | `--dat` bulk flat-file (or `--query` REST) |
| `uniref` | UniRef{50,90,100} deduped protein clusters | `--file` bulk FASTA, `--stride` to subsample |
| `ensembl` | Ensembl peptides (rich FASTA headers) | bulk FASTA |
| `ensembl-gff` | gene models (biotype, location, exon/transcript counts) | bulk GFF3 |
| `ensembl-regulatory` | Regulatory-Build features (enhancer/promoter/CTCF/…) + DNA | `--genome` local (or REST) |
| `ensembl-splice` | 5′ donor / 3′ acceptor junction windows + DNA | `--genome` local (or REST) |
| `ensembl-dogma` | **DNA (pre-mRNA) + spliced RNA + protein** for one transcript | `--genome` + FASTAs local (or REST) |

Ensembl sources are **multi-species**: `--species` is a name, comma-list, or
`all` (~356 vertebrates), resolved to per-species URLs + organism + taxid.

## Design (holds for every record)

- **Sequence-first.** `>id name [organism]` → `<seq>…</seq>` → metadata. In an
  autoregressive LM this trains P(metadata | sequence) — recognition, the useful
  inference direction. `--ordering both` also emits the design direction.
- **Provenance.** `source, source_version, source_url, license, accession,
  organism, taxid` on every record (TheBioCollection shipped only `{text, record_type}`).
- **Verified, not asserted.** Dogma checks `translate(CDS)==protein` and
  `mRNA==spliced exons`; splice reports/classifies the *observed* motif; records
  that fail are skipped, not rendered wrong.
- **Clean & deduped.** Missing fields are omitted (no empty-slot bugs);
  exact-sequence dedup within a run.
- **Bulk vs REST.** Billions of tokens come from local bulk files; REST is for
  one-species pilots.

## Record schema (one JSON object per line)

`id · source · source_version · source_url · license · accession · entity_type`
(`protein|transcript|gene|regulatory_feature|splice_junction|central_dogma`) ·
`seq_type` (`aa|dna|rna`) · `seq_len · organism · taxid · gene · name ·
annotations{…} · sequence · sequences{dna_genomic,rna,cds,protein}` (dogma only)
`· ordering · text`. `text` is what a tokenizer consumes; the rest is queryable.

## Example records

The `text` field, verbatim from [`samples/`](samples/) — sequences and long lists
truncated with `…`.

**Protein** (Swiss-Prot; UniProt annotation rendered as prose):

```
>sp:P69905 Hemoglobin subunit alpha [Homo sapiens (Human)]
<protein>MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPH…</protein>

Hemoglobin subunit alpha — UniProtKB/Swiss-Prot P69905 — is a 142-residue protein from Homo sapiens (Human) (NCBI taxon 9606); gene HBA1.
Function: Involved in oxygen transport from the lung to the various peripheral tissues
GO annotations: blood microparticle (component); hemoglobin complex (component); heme binding (function); oxygen binding (function); …
Sequence features: 151 variant; 17 mod res [4 (Phosphoserine)]; 12 helix; 2 binding [59; 88 (proximal binding residue)]; …
Keywords: 3D-structure; Heme; Hereditary hemolytic anemia; Iron; Metal-binding; Oxygen transport; …
Cross-references: PDB 1A00, 1A01, 1A0U +343 more; Reactome R-HSA-1237044 …; Pfam PF00042
Lineage: Eukaryota > Metazoa > Chordata > … > Mammalia > Eutheria
```

**UniRef50 cluster** (the deduped protein backbone — compact by design):

```
>uniref50:A0A007 MoeK5 [Streptomyces viridosporus]
<protein>MGYIHTALKSAGFHHVIQVDTPALGLDSEGLRKLLADFEPDLVGVS…</protein>

MoeK5 — UniRef50 A0A007 — is a 407-residue protein from Streptomyces viridosporus (NCBI taxon 67581).
Cluster size: 1 member sequence(s)
Representative: A0A007_STRVD
```

**Central dogma** — DNA → RNA → protein in one document, verified:

```
>ensembl-dogma:ENST00000641515 OR4F5 [Homo sapiens]
Gene OR4F5 (protein_coding), Ensembl canonical transcript ENST00000641515, 1:65,419-71,585 (+), 3 exon(s) / 2 intron(s).

Genomic DNA (pre-mRNA, sense strand; exons UPPERCASE, introns lowercase):
<dna>CCCAGATCTCTTCAGgtacatctagtccattcataaagggcttttaattaaccaag…</dna>

Transcription and splicing remove 2 intron(s) (3549 nt) to give the mature mRNA (2618 nt):
<rna>CCCAGAUCUCUUCAGUUUUUAUGCCUCAUUCUGUGAAAAUUGCUGUAGUCUCUUCC…</rna>
Exon boundaries in the mRNA: 1-15, 16-69, 70-2618.
5' UTR: 1-60 (60 nt); CDS: 61-1041 (981 nt); 3' UTR: 1042-2618 (1577 nt).

Translation of the CDS (981 nt, standard genetic code) yields the 326-residue protein:
<protein>MKKVTAEAISWNESTSETNNSMVTEFIFLGLSDSQELQTFLFMLFF…</protein>

Verified: the mRNA equals the genomic exons with introns removed; translate(CDS) equals the protein.
```

The splice is legible in the DNA itself: `…TTCAG` (exon, UPPER) → `gtacatct…`
(intron, lower, canonical `gt`).

**Splice site** (donor; the acceptor record is its mirror, intron→exon):

```
>ensembl:ENST00000456328:intron1:donor [Homo sapiens]
<dna>CCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCATAGGGG…</dna>

The sequence is a 200 bp window over a 5' splice donor site in gene DDX11L2 from Homo sapiens (NCBI taxon 9606).
Splice site: 5' splice donor
Transcript: ENST00000456328
Intron: 1 of 2
Intron location: 1:12,228-12,612 (+)
Intron length: 385 bp
Donor dinucleotide: GT — canonical (GT-AG)
```

Gene-model, regulatory-feature, multi-species and TrEMBL examples are in
[`samples/`](samples/).

## Usage

```bash
# proteins (bulk, local)
python build_bio_corpus.py uniprot --dat uniprot_sprot.dat.gz --ordering both --out sprot.jsonl
python build_bio_corpus.py uniref  --file uniref50.fasta.gz --stride 3 --out uniref50.jsonl

# central dogma — REST for a quick species sample, or --genome for offline scale
python build_bio_corpus.py ensembl-dogma --view all --min-exons 2 --limit 20 --out dogma.jsonl
python build_bio_corpus.py ensembl-dogma --species mus_musculus --gff mouse.gff3.gz \
    --genome mouse.dna.toplevel.fa.gz --cdna mouse.cdna.fa.gz --cds mouse.cds.fa.gz --pep mouse.pep.fa.gz \
    --view all --min-exons 2 --limit 0 --out dogma_mouse.jsonl

# splice / regulatory (add --genome for offline)
python build_bio_corpus.py ensembl-splice --site both --limit 20 --out splice.jsonl
python build_bio_corpus.py ensembl-regulatory --limit 20 --out reg.jsonl

# multi-species proteins
python build_bio_corpus.py ensembl --species "mus_musculus,danio_rerio,gallus_gallus" --limit 50 --out fish_etc.jsonl

# de-duplicate across shards (priority order: Swiss-Prot's copy of a protein beats UniRef's beats Ensembl's)
python dedup.py sprot.jsonl uniref50.jsonl ensembl_pep.jsonl dogma_*.jsonl --out corpus.jsonl
```

`--limit` is per species (default 50; `--limit 0` = uncapped).

## Scale & token budget

Under Marin's tokenizer (Llama-3): protein ≈ **0.55 tok/residue**. Rough totals —
Swiss-Prot ~1 B (both orderings), **UniRef50 ~12.6 B full** (`--stride` to dial),
human genomic/dogma/splice/regulatory ~0.3 B. Assemble a **4–10 B corpus** with
**[RECIPE.md](RECIPE.md)**; full method + per-subset numbers in
**[TOKEN_ESTIMATE.md](TOKEN_ESTIMATE.md)** (reproduce with `estimate_tokens.py`).

## Samples & caveats

- [`samples/`](samples/): committed example output for every source.
- **Dedup**: exact-sequence within a run; `dedup.py` (priority order) is the
  cross-file pass that collapses the same protein across Swiss-Prot / UniRef /
  Ensembl (namespaced by seq type, so dogma records are never dropped).
- **Regulatory** features exist for human/mouse only (Ensembl builds them there).
- **REST mode** is bounded to one species / a pilot; use `--genome` for scale.
- **Nucleotide is token-expensive** (~0.55 tok/nt) — prefer gene/feature records
  over raw genome dumps; cap with `--max-dna`. This is a POC, not a pipeline
  (no distributed execution / cross-run dedup / bulk TrEMBL).
