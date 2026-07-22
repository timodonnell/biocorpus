# Corpus recipe — targeting 4–10 B tokens

At this scale the **REST/pagination path doesn't work** (billions of records =
millions of requests). The 4–10 B has to come from **bulk files streamed
locally**. The tunable bulk is **UniRef50** (deduped proteins); everything else
is a fixed accent.

Token rates below are measured with Marin's tokenizer (Llama-3) — see
[TOKEN_ESTIMATE.md](TOKEN_ESTIMATE.md).

## Two regimes

| regime | sources | scales to billions? |
|---|---|---|
| **bulk file, local** | `uniprot --dat`, `uniref --file`, `ensembl-* --genome` | **yes** — download once, stream/seek locally |
| REST-backed | `ensembl-*` (default), `uniprot --query` | bounded — fine for one species / a pilot, not billions |

## Recommended ~6 B recipe

| # | component | command | records | tok/rec | **tokens** |
|--:|---|---|--:|--:|--:|
| 1 | **Swiss-Prot** (both orderings) | `uniprot --dat uniprot_sprot.dat.gz --ordering both --limit 0` | 1.15 M | ~900 | **~1.0 B** |
| 2 | **UniRef50** (stride 3) | `uniref --file uniref50.fasta.gz --stride 3 --limit 0` | 12.9 M | ~325 | **~4.2 B** |
| 3 | Human genomic (dogma+splice) | `ensembl-dogma … ` / `ensembl-splice …` (human) | ~0.4 M | ~1–2 k | **~0.5 B** |
| 4 | Human regulatory | `ensembl-regulatory --limit 0` | 612 K | ~330 | **~0.2 B** |
| | | | | **total** | **≈ 5.9 B** |

## Tuning: UniRef50 `--stride` is the dial

Full UniRef50 (38.8 M clusters) is **~12.6 B** on its own, so subsample it:

| target total | UniRef50 stride | UniRef50 clusters | UniRef50 tokens | + Swiss-Prot + genomic |
|---|--:|--:|--:|---|
| **~4 B** | 5 | 7.8 M | ~2.5 B | +1.0 +0.5 = **~4.0 B** |
| **~6 B** | 3 | 12.9 M | ~4.2 B | +1.0 +0.7 = **~5.9 B** |
| **~8 B** | 2 | 19.4 M | ~6.3 B | +1.0 +0.7 = **~8.0 B** |
| **~10 B** | 2 + UniRef90 top-up, or full UniRef50 minus genomic | — | — | **~10 B** |

`--stride N` emits every Nth cluster; because the UniRef FASTA is sorted by
descending length, striding gives a **length-stratified, representative** subset
(plain `--limit` alone would grab only the giant-protein head — 40 k+ aa each).

## Full command sequence

```bash
BASE=https://ftp.uniprot.org/pub/databases/uniprot
DL () { curl -sL "$1" -o "$2"; }

# --- bulk downloads (one-time) ---
DL $BASE/current_release/knowledgebase/complete/uniprot_sprot.dat.gz uniprot_sprot.dat.gz   # ~150 MB
DL $BASE/uniref/uniref50/uniref50.fasta.gz                           uniref50.fasta.gz      # 8.8 GB

# --- 1. Swiss-Prot, both orderings (~1.0 B) ---
python build_bio_corpus.py uniprot --dat uniprot_sprot.dat.gz --ordering both --limit 0 --out out/sprot.jsonl

# --- 2. UniRef50 deduped backbone, tuned to the token budget (~4.2 B at stride 3) ---
python build_bio_corpus.py uniref --file uniref50.fasta.gz --stride 3 --limit 0 --out out/uniref50.jsonl

# --- 3–4. Human genomic accent (REST; runs for hours — background it) ---
python build_bio_corpus.py ensembl-dogma      --view all --min-exons 2 --limit 0 --out out/dogma_human.jsonl
python build_bio_corpus.py ensembl-splice     --site both             --limit 0 --out out/splice_human.jsonl
python build_bio_corpus.py ensembl-regulatory                         --limit 0 --out out/reg_human.jsonl

# optional: a few model vertebrates (multiplies REST time)
python build_bio_corpus.py ensembl-dogma --species "mus_musculus,danio_rerio,gallus_gallus,xenopus_tropicalis" \
    --view all --min-exons 2 --limit 0 --out out/dogma_species.jsonl
```

`--limit 0` = uncapped (the default `--limit` is 50, for quick sampling).

## Run notes & honest limits

- **Steps 1–2 are the corpus** (~5 B) and are pure local streaming: fast, no REST,
  reproducible. Output is ~15–20 GB of JSONL; pipe through `zstd` if you want
  `.jsonl.zst` shards like TheBioCollection.
- **Genomic (steps 3–4) scales offline.** Default is REST (fine for one species /
  a pilot). At scale pass `--genome <species>.dna.toplevel.fa[.gz]` to read DNA
  from the **local genome FASTA** (indexed seek; validated byte-identical to REST)
  — no per-record requests. Dogma also takes `--cdna/--cds/--pep` (default:
  streamed once per species). Loop species in a shell script, downloading each
  genome, to go all-vertebrates fully offline:

  ```bash
  SP=mus_musculus ASM=GRCm39 REL=112; Sp=$(python -c "print('$SP'.capitalize())")
  base=https://ftp.ensembl.org/pub/release-$REL
  for k in dna/$Sp.$ASM.dna.toplevel cdna/$Sp.$ASM.cdna.all cds/$Sp.$ASM.cds.all pep/$Sp.$ASM.pep.all; do
    curl -sL $base/fasta/$SP/$k.fa.gz -o $(basename $k).fa.gz; done
  curl -sL $base/gff3/$SP/$Sp.$ASM.$REL.gff3.gz -o gff.gff3.gz
  python build_bio_corpus.py ensembl-dogma --species $SP --gff gff.gff3.gz \
    --genome $Sp.$ASM.dna.toplevel.fa.gz --cdna $Sp.$ASM.cdna.all.fa.gz \
    --cds $Sp.$ASM.cds.all.fa.gz --pep $Sp.$ASM.pep.all.fa.gz \
    --view all --min-exons 2 --limit 0 --out out/dogma_$SP.jsonl
  python build_bio_corpus.py ensembl-splice --species $SP --gff gff.gff3.gz \
    --genome $Sp.$ASM.dna.toplevel.fa.gz --site both --limit 0 --out out/splice_$SP.jsonl
  ```
- **Redundancy:** UniRef50 already removes near-duplicate proteins; if you add
  multi-species Ensembl proteins/dogma, expect ortholog overlap with UniRef50 —
  dedup by sequence (the builder does exact-sequence dedup within a run; a
  cross-run pass would catch the rest).
- **Both orderings** doubles document count for the components you apply it to
  (Swiss-Prot here); it does not change which sequences are covered.
