# Sequence-first biology pre-training corpus — proof of concept

A small, runnable builder that turns **UniProtKB/Swiss-Prot** and **Ensembl**
entries into pre-training documents where the **sequence comes first** and the
metadata follows. It exists to make one design idea concrete and testable, and
to fix — by construction — the specific problems found in the
[TheBioCollection review](../../analyses/thebiocollection-review.html).

```
>sp:P69905 Hemoglobin subunit alpha [Homo sapiens (Human)]
<protein>MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG…</protein>

Hemoglobin subunit alpha — UniProtKB/Swiss-Prot P69905 — is a 142-residue protein from Homo sapiens (Human) (NCBI taxon 9606); gene HBA1.
Function: Involved in oxygen transport from the lung to the various peripheral tissues
Involvement in disease: Heinz body anemias (HEIBAN) [MIM:140700]: Form of non-spherocytic hemolytic anemia of Dacie type 1. …
GO annotations: blood microparticle (component); heme binding (function); oxygen binding (function); …
Sequence features: 151 variant; 17 mod res [4 (Phosphoserine); …]; 12 helix; 2 binding [59; 88 (proximal binding residue)]; …
Keywords: 3D-structure; Heme; Hereditary hemolytic anemia; Iron; Metal-binding; Oxygen transport; …
Cross-references: PDB 1A00, 1A01, 1A0U, 1A0Z +343 more; Reactome R-HSA-1237044, …; Pfam PF00042
Lineage: Eukaryota > Metazoa > Chordata > … > Mammalia > Eutheria
```

## Why sequence-first

An autoregressive LM factorises its loss left-to-right, so **the order you write
a record decides which conditional the model is trained on.**

| Ordering | Trains | Capability |
|---|---|---|
| `sequence_first` (default) | P(metadata \| sequence) | **recognition** — "given a sequence, tell me what it is" |
| `metadata_first` | P(sequence \| metadata) | **design** — "given a description, generate a sequence" |

For a science assistant, recognition is the direction you actually query at
inference, so sequence-first is the right default. (TheBioCollection's dominant
73% molecule records are already SMILES-first — recognition-ordered — and
molecule-property recognition was one of its real eval gains. This extends that
pattern to proteins and genes.) Use `--ordering both` to emit one document of
each and get both directions.

A short identity **anchor** (`>db:accession name [organism]`) precedes the
sequence so the sequence is never fully context-free, which helps the model bind
sequence ↔ identity.

## What it fixes from TheBioCollection

| Review finding | Fix here |
|---|---|
| **#10 No provenance** — records were only `{text, record_type}` | Every record carries `source`, `source_version`, `source_url`, `license`, `accession`, `organism`, `taxid`, plus the structured `annotations`. |
| **#4 Empty-slot template bugs** (~5% of records: `performed at pH  and .`) | Missing fields are **omitted**, never rendered blank. |
| **#8 Redundancy** | Exact-sequence dedup (a stand-in for UniRef clustering); records are dropped if a sequence is already seen. |
| **#9 Loose semantics** | `entity_type`, `seq_type`, and typed `annotations` instead of a filename-echoing `record_type`. |
| **#2 Thin, templated prose** | Metadata rendered as clean natural language from curated fields (function, disease, GO, features, x-refs), not raw flat-files. |

## Record schema

One JSON object per line. `text` is what a tokenizer consumes; everything else
is queryable metadata.

| Field | Type | Notes |
|---|---|---|
| `id` | str | globally unique, e.g. `uniprot:P69905` (suffixed `#se`/`#me` under `--ordering both`) |
| `source` | str | `uniprot_swissprot` \| `ensembl` |
| `source_version` | str | release label (e.g. `release-112`, `current_release`) |
| `source_url` | str | exact URL the entry came from |
| `license` | str | UniProt = CC-BY-4.0; Ensembl = open (verify assembly terms) |
| `accession` | str | primary accession / stable id |
| `entity_type` | str | `protein` \| `transcript` |
| `seq_type` | str | `aa` \| `dna` |
| `seq_len` | int | length in residues/nt |
| `organism`, `taxid` | str | species + NCBI taxon id |
| `gene`, `name` | str | gene symbol; display name |
| `annotations` | object | source-specific typed fields (function, go, features, keywords, xrefs, lineage, description, location, biotype, …) |
| `sequence` | str | raw residues/bases |
| `ordering` | str | `sequence_first` \| `metadata_first` |
| `text` | str | the rendered pre-training document |

## Usage

```bash
# Rich, curated protein records (fetched individually via the UniProt REST API)
python build_bio_corpus.py uniprot \
  --accessions P02768,P00533,P69905,P0DTC2,P00698 \
  --ordering sequence_first --out samples/uniprot_swissprot.sample.jsonl

# Stream the head of the full Swiss-Prot flat file instead (no full download)
python build_bio_corpus.py uniprot --limit 200 --out out.jsonl

# Ensembl peptides (default organism: yeast). Use --seq-type cdna for nucleotide transcripts.
python build_bio_corpus.py ensembl --limit 40 \
  --ordering sequence_first --out samples/ensembl_pep.sample.jsonl

# Ensembl GFF3 gene models, joined to the canonical protein product
# (--seq-type none emits annotation-only gene records)
python build_bio_corpus.py ensembl-gff --limit 40 \
  --ordering sequence_first --out samples/ensembl_gff.sample.jsonl

# Sample a representative slice of Swiss-Prot or TrEMBL via the UniProt query API
python build_bio_corpus.py uniprot --query "reviewed:true"  --limit 500 --out sprot.jsonl
python build_bio_corpus.py uniprot --query "reviewed:false" --limit 500 --out trembl.jsonl

# Emit both recognition and design orderings
python build_bio_corpus.py uniprot --accessions P69905 --ordering both --out both.jsonl
```

The `ensembl-gff` source produces **gene-level records** (`entity_type: gene`)
carrying the GFF-derived gene model — biotype, genomic location, and
transcript/exon counts — joined by `gene_id` to the canonical product sequence:

```
>ensembl:YAL067C SEO1 [Saccharomyces cerevisiae]
<protein>MYSIVKEIIVDPYKRLKWGFIPVKRQVEDLPDDLNSTEIVTISNSIQSHETAENF…</protein>

SEO1 — Ensembl gene YAL067C — is a protein_coding gene from Saccharomyces cerevisiae (NCBI taxon 4932).
Gene model: 1 transcript(s), 1 exon(s) in the canonical transcript
Canonical product: canonical protein product, 593 residues
Description: Putative permease; member of the allantoate transporter subfamily of the major facilitator superfamily; …
Location: I:7,235-9,016 (-)
Biotype: protein_coding
```

Requires `biopython` (Swiss-Prot parsing) and stdlib only otherwise. Committed
example output is in [`samples/`](samples/): `uniprot_swissprot.sample.jsonl`
(proteins), `ensembl_pep.sample.jsonl` (peptides from FASTA headers), and
`ensembl_gff.sample.jsonl` (gene models from GFF3 with the canonical product
sequence and transcript/exon counts).

## Caveats — this is a POC, not a pipeline

A production build would add:

- **Real dedup/clustering:** UniRef50/90 for proteins instead of exact-match; genome-window dedup for nucleotides.
- **Richer Ensembl records:** GFF3 gene models (biotype, transcript/exon counts, location) are built and joined to the canonical product; a production build would add transcript-level cDNA records, GO terms and orthologs (BioMart/REST), and would still **avoid raw genome dumps** — nucleotide sequence is long and low-information per token, so prefer gene/transcript-level records and cap `--max-seq-len`.
- **Tokenization:** if the bio fraction grows, add sequence-aware or byte-level tokenization; raw AA/nt in a text BPE vocab is inefficient.
- **Eval decontamination:** dedup against downstream biology benchmarks (including TheBioCollection-Eval) before training.
- **Scale + execution:** streaming/sharding, resumability, and packing many short records per document (cf. Marin's `bio_chem` datakit, which does exactly this for capped pilot slices).
- **Mixture discipline:** keep bio a small single-digit % of a general corpus unless training a bio-specialised model.

## Sources & licensing

- **UniProtKB/Swiss-Prot** — [CC-BY 4.0](https://www.uniprot.org/help/license). Streamed from the UniProt FTP / REST API.
- **Ensembl** — annotation is distributed with [no restrictions](https://www.ensembl.org/info/about/legal/disclaimer.html); some genome assemblies carry third-party terms, so verify per species.
