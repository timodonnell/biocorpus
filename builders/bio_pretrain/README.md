# Sequence-first biology pre-training corpus builder

Turns biological databases into **sequence-first** pre-training documents — raw
sequence, then clean natural-language annotation — with provenance, verified
mappings, de-duplication, and an **annotation-quality gate**. Built to fix the
failure modes found in the
[TheBioCollection review](../../analyses/thebiocollection-review.html): no
provenance, ~5% empty-slot template bugs, "verifiable" tasks that aren't
learnable, and DNA/RNA/protein never co-occurring in one document.

`build_bio_corpus.py` renders records; `uniref_filter.py` + `run_v2.sh` assemble
the full corpus; `dedup.py` collapses duplicates. Only dependency is `biopython`.

The current build (**biocorpus v2**) is **21.1 M documents / ~9.15 B tokens**
(Llama-3 / `marin-community/marin-tokenizer`), published at
[hf://buckets/timodonnell/biocorpus](https://huggingface.co/buckets/timodonnell/biocorpus).

| source | documents | tokens |
|---|--:|--:|
| TrEMBL representatives (annotated, UniRef50-deduped) | 20,091,178 | 8.49 B |
| Swiss-Prot (both orderings) | 950,174 | 0.58 B |
| human central-dogma + splice | 96,768 | 0.07 B |

## Sources

| subcommand | what | scale path (billions) |
|---|---|---|
| `uniprot` | Swiss-Prot / TrEMBL proteins (function, GO, features) | `--dat` bulk flat-file (shards from `uniref_filter.py`) |
| `uniref` | UniRef{50,90,100} deduped protein clusters | `--file` bulk FASTA, `--stride` to subsample |
| `ensembl` | Ensembl peptides (rich FASTA headers) | bulk FASTA |
| `ensembl-gff` | gene models (biotype, location, exon/transcript counts) | bulk GFF3 |
| `ensembl-regulatory` | Regulatory-Build features (enhancer/promoter/CTCF/…) + DNA | `--genome` local (or REST) |
| `ensembl-splice` | 5′ donor / 3′ acceptor junction windows + DNA | `--genome` local (or REST) |
| `ensembl-dogma` | **DNA (pre-mRNA) + spliced RNA + protein** for one transcript | `--genome` + FASTAs local (or REST) |

Ensembl sources are **multi-species**: `--species` is a name, comma-list, or
`all` (~356 vertebrates), resolved to per-species URLs + organism + taxid.

## The protein backbone: annotated + deduplicated (no REST)

The bulk of the corpus is the **UniProtKB entry of each UniRef50 cluster
representative** — the deduplicated protein universe (50% identity, ~34 M
UniProtKB representatives), each carrying its real annotation. Obtained by a
local flat-file join, no per-record web requests:

```bash
# phase 1 — subset the TrEMBL flat file to UniRef50 representatives (fast byte-level filter)
python uniref_filter.py --in uniprot_trembl.dat.gz --uniref uniref50.fasta.gz \
    --repset-cache repset50.pkl --out-prefix out/trembl_rep --shards 32
# phase 2 — render each shard with the ordinary uniprot source (gate + dedup apply)
python build_bio_corpus.py uniprot --dat out/trembl_rep.00.dat --ordering sequence_first --limit 0 --out out/trembl_rep.00.jsonl
```

`run_v2.sh` wires the whole build end-to-end (Swiss-Prot + TrEMBL reps + human
genomic → dedup). See [RECIPE.md](RECIPE.md) for download URLs and knobs.

## Annotation-quality gate

A protein is kept **only** if it carries genuinely learnable annotation:

- **tier A (biological):** function, catalytic activity, subcellular location, pathway, disease, or GO; **or**
- **tier B (structural):** a domain boundary, active/binding site, signal peptide, transmembrane span, PTM, disulfide, repeat, …

Entries whose only annotation is a disordered/coiled-coil region — most
"Uncharacterized protein" TrEMBL entries — are dropped (`--allow-unannotated`
disables the gate). Genomic records are intrinsically annotated and always kept.
UniProt evidence tags (`{ECO:…}`), empty `Evidence=` clauses, and markup are
stripped; no `<protein>`/`<dna>` tags.

## Design (holds for every record)

- **Sequence-first.** `>id name [organism]` → sequence → annotation. In an
  autoregressive LM this trains P(annotation | sequence) — recognition, the useful
  inference direction. `--ordering both` also emits the design direction (Swiss-Prot).
- **Provenance.** `source, source_version, source_url, license, accession,
  organism, taxid` on every record (TheBioCollection shipped only `{text, record_type}`).
- **Verified, not asserted.** Dogma checks `translate(CDS)==protein` and
  `mRNA==spliced exons`; splice reports the *observed* motif; records that fail are skipped.
- **Clean & deduped.** Missing fields omitted (no empty-slot bugs); exact-sequence
  dedup within a run and across files (`dedup.py`).
- **Bulk vs REST.** Billions of tokens come from local bulk files; REST is for one-species pilots.

## Record schema (one JSON object per line)

`id · source · source_version · source_url · license · accession · entity_type`
(`protein|transcript|gene|regulatory_feature|splice_junction|central_dogma`) ·
`seq_type` (`aa|dna|rna`) · `seq_len · organism · taxid · gene · name ·
annotations{…} · sequence · sequences{dna_genomic,rna,cds,protein}` (dogma only)
`· ordering · text`. `text` is what a tokenizer consumes; the rest is queryable.

## Example records (v2 `text`, verbatim)

**TrEMBL representative** (tier A — deduplicated, richly annotated):

```
>tr:Q977Q6 Methionine aminopeptidase [uncultured crenarchaeote 4B7]
MTFDNYIKAGKIAGEIRENVRKTDWVGKTVYEICEYVENEIKKRGAKCAFPVNTSINEVAAHYTAEPNDEIT…

Methionine aminopeptidase — UniProtKB/TrEMBL Q977Q6 — is a 225-residue protein from uncultured crenarchaeote 4B7 (NCBI taxon 44557).
Catalytic activity: Reaction=Release of N-terminal amino acids, preferentially methionine…; EC=3.4.11.18
GO annotations: cytoplasm (component); initiator methionyl aminopeptidase activity (function); metal ion binding (function); metalloexopeptidase activity (function); proteolysis (process)
Sequence features: 1 domain [6-194 (Peptidase M24)]
Keywords: Aminopeptidase; Hydrolase; Metal-binding; Protease
Lineage: Archaea > Nitrososphaerota > Nitrososphaeria > Nitrosopumilales > environmental samples
```

**TrEMBL representative** (tier B — no biological field, kept on a domain boundary):

```
>tr:A0ABD5XRM9 Bacterio-opsin activator domain-containing protein [Halobaculum litoreum]
MGVHAAVTVRAREFALARTLAVAPSARVTLEPVVPFGAGFAPAVRIRADDPDLVTDLVAAEADVRAVEP…

Bacterio-opsin activator domain-containing protein — UniProtKB/TrEMBL A0ABD5XRM9 — is a 178-residue protein from Halobaculum litoreum (NCBI taxon 3031998); gene ACFQRB_18020.
Sequence features: 1 domain [7-133 (Bacterioopsin transcriptional activator GAF and HTH associated)]; 1 region [132-178 (Disordered)]
Keywords: Reference proteome
Lineage: Archaea > Methanobacteriati > Methanobacteriota > Stenosarchaea group > Halobacteria > Halobacteriales > Haloferacaceae > Halobaculum
```

**Central dogma** — DNA → spliced RNA → protein in one document, verified:

```
>ensembl-dogma:ENST00000633214 OVOL3 [Homo sapiens]
Gene OVOL3 (protein_coding), Ensembl canonical transcript ENST00000633214, 19:36,111,143-36,113,711 (+), 4 exon(s) / 3 intron(s).

Genomic DNA (pre-mRNA, sense strand; exons UPPERCASE, introns lowercase):
GGGCTGAGGTCTGACAGCAGGTGGAAGCAGCCCCTGTGTGTGGAGAGCCTTCCGGAGGGCATGCCCCGCGCC…CCCAGgtgggcccctcactgtgcctggagg…

Transcription and splicing remove 3 intron(s) (1886 nt) to give the mature mRNA (683 nt):
GGGCUGAGGUCUGACAGCAGGUGGAAGCAGCCCCUGUGUGUGGAGAGCCUUCCGGAGGG…
Exon boundaries in the mRNA: 1-154, 155-219, 220-424, 425-683.
5' UTR: 1-60 (60 nt); CDS: 61-633 (573 nt); 3' UTR: 634-683 (50 nt).

Translation of the CDS (573 nt, standard genetic code) yields the 190-residue protein:
MPRAFLVRSRRPQPPNWGHLPDQLRGDAYIPDCSSLGGPPAQQSSSVRDPWTAQPTQGNLTSAP…

Verified: the mRNA equals the genomic exons with introns removed; translate(CDS) equals the protein.
```

**Splice site** (donor; the acceptor record is its mirror, intron→exon):

```
>ensembl:ENST00000594440:intron1:donor [Homo sapiens]
CCTTGATTAAACGTGCACTTCGCAGTCCTCGGTTCTCCATACCCGTGACCTGGGGATCGCTACGGACCTT…

The sequence is a 200 bp window over a 5' splice donor site in gene ZNF841 from Homo sapiens (NCBI taxon 9606).
Splice site: 5' splice donor
Transcript: ENST00000594440
Intron: 1 of 6
Intron length: 1,603 bp
Donor dinucleotide: GT — canonical (GT-AG)
```

Gene-model, regulatory-feature and Swiss-Prot examples are in [`samples/`](samples/).

## Usage

```bash
# full v2 build (downloads listed in RECIPE.md must be present)
bash run_v2.sh

# individual sources
python build_bio_corpus.py uniprot --dat uniprot_sprot.dat.gz --ordering both --limit 0 --out sprot.jsonl
python build_bio_corpus.py ensembl-dogma --species mus_musculus --gff mouse.gff3.gz \
    --genome mouse.dna.toplevel.fa.gz --cdna mouse.cdna.fa.gz --cds mouse.cds.fa.gz --pep mouse.pep.fa.gz \
    --view all --min-exons 2 --limit 0 --out dogma_mouse.jsonl

# cross-file dedup (priority order: Swiss-Prot beats TrEMBL reps beats genomic)
python dedup.py sprot.jsonl trembl_rep.*.jsonl dogma_*.jsonl splice_*.jsonl --out corpus.jsonl
```

`--limit` is per species (default 50; `--limit 0` = uncapped).

## Scale & token budget

Under Marin's tokenizer (Llama-3): protein ≈ **0.55 tok/residue**. The v2 corpus
is **9.15 B tokens** (annotated TrEMBL reps 8.49 B + Swiss-Prot 0.58 B + human
genomic 0.07 B). Dial the protein bulk with the UniRef identity level and the
quality gate; assemble variants with **[RECIPE.md](RECIPE.md)** (per-subset numbers
in **[TOKEN_ESTIMATE.md](TOKEN_ESTIMATE.md)**, reproduce with `estimate_tokens.py`).

## Caveats

- **Regulatory** features exist for human/mouse only (Ensembl builds them there).
- **REST mode** is bounded to one species / a pilot; use `--genome`/`--dat` for scale.
- **Nucleotide is token-expensive** (~0.55 tok/nt) — prefer gene/feature records
  over raw genome dumps; cap with `--max-dna`.
- Genomic slice is currently human only; the join is single-node (no distributed
  execution). The annotated-UniRef path uses a local TrEMBL flat file (~118 GB gz).
