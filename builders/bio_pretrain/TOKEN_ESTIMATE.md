# Token-yield estimate (Marin tokenizer)

How many tokens would a sequence-first UniProt/Ensembl corpus produce, per
database subset, under **Marin's tokenizer** — `marin-community/marin-tokenizer`,
which is the **Llama-3 128,256-vocab** tokenizer (confirmed in Marin's
`experiments/marin_tokenizer.py`).

Reproduce with [`estimate_tokens.py`](estimate_tokens.py):

```bash
hf download marin-community/marin-tokenizer tokenizer.json --local-dir ./marin_tok
python build_bio_corpus.py uniprot --query "reviewed:true"  --limit 1000 --out sprot.jsonl
python build_bio_corpus.py uniprot --query "reviewed:false" --limit 1000 --out trembl.jsonl
python estimate_tokens.py --tokenizer ./marin_tok/tokenizer.json --db-mean-aa 330 \
    --sample "Swiss-Prot=sprot.jsonl" --sample "TrEMBL=trembl.jsonl" \
    --sample "Ensembl(pep)=samples/ensembl_pep.sample.jsonl"
```

## Tokenization efficiency (the key robust number)

| Content | chars/token | tokens per unit |
|---|--:|--:|
| Protein sequence (amino acids) | 1.81 | **0.55 tok/residue** |
| Nucleotide sequence (DNA/cDNA) | 1.76 | 0.57 tok/nt |
| English prose (metadata) | 4.29 | — |

Sequences are **~2.4× less token-efficient per character** than prose: a 330-aa
protein ≈ 180 tokens of low-information sequence. This is stable across every
sample (0.55 ± 0.005 tok/residue) and is organism-independent, so it anchors the
projection below.

## Measured per-record (samples of ~1,000, rendered sequence-first)

| sample | n | mean aa | seq tok | meta tok | total tok |
|---|--:|--:|--:|--:|--:|
| Swiss-Prot (reviewed) | 1000 | 555 | 306 | 732 | 1038 |
| TrEMBL (unreviewed) | 1000 | 779 | 429 | 586 | 1014 |
| Ensembl human (pep) | 30 | 351 | 194 | 140 | 333 |

⚠️ The UniProt search API returns well-annotated (and longer) entries first, so
these sample **lengths and metadata run high** vs. the true database means
(Swiss-Prot ~355 aa; TrEMBL ~330 aa with many fragments). The projection below
therefore anchors the sequence term on a realistic 330-aa mean × 0.55, and
treats metadata as a band.

## Projected whole-database tokens

Entry counts fetched 2026-07 from the UniProt / UniRef REST APIs.

| Subset | entries | tok/record | **total tokens** | basis |
|---|--:|--:|--:|---|
| **Swiss-Prot** (reviewed) | 575,503 | ~1000 | **~0.6 B** | measured (rich) |
| **TrEMBL** (unreviewed) | 149,234,636 | 207 → 1014 | **~31 B (floor) – 151 B (ceiling)** | seq-floor → rich sample |
| UniRef100 | 220,919,788 | 207 | ~46 B | seq + identity line |
| UniRef90 | 121,389,642 | 207 | ~25 B | seq + identity line |
| **UniRef50** | 38,794,121 | 207 | **~8 B** | seq + identity line |
| Ensembl human (protein-coding genes) | ~20,000 | 333 | ~7 M | measured (pep) |

`tok/record = 207` is `330 aa × 0.55 tok/residue + ~25` for a minimal identity
line; realistic TrEMBL metadata pushes it toward the middle of the 31–151 B band
(most likely **~40–70 B**). Emitting **both** orderings (recognition + design)
doubles every figure.

## What this means for the mixture

For scale reference, Marin 8B trained on **>12 T** tokens.

- **Swiss-Prot (~0.6 B) is the obvious first ingredient** — small, uniformly
  high-quality, richly annotated, permissive (CC-BY). At <0.01% of a 12 T budget
  it is essentially free; include it (and probably at both orderings).
- **Raw TrEMBL (~31–151 B) is not worth it as-is** — it's dominated by
  low-information, highly redundant sequence tokens. If you want protein-sequence
  breadth, **dedup to UniRef50 (~8 B)** — a ~4–18× reduction that removes near-duplicate
  homologs while keeping diversity. UniRef50 is the natural "protein backbone".
- **Ensembl protein content overlaps UniProt** and one species is negligible
  (~7 M). Ensembl's real value is **genomic/regulatory** (gene models, splice
  sites, regulatory regions) — but nucleotide sequence is even less
  token-efficient, so prefer **gene/annotation records over raw genome dumps**
  and cap sequence length.
- **Sequences are expensive tokens.** ~50–90% of every record here is sequence
  tokens carrying little signal for a *text* model. That is the strongest argument
  for (a) keeping the bio fraction small, (b) deduping hard (UniRef50), and
  (c) considering sequence-aware tokenization if the fraction grows.

**Rule of thumb:** Swiss-Prot (~0.6 B, both orderings ≈ 1.2 B) + UniRef50 (~8 B)
+ curated Ensembl gene/annotation records is a ~10 B-token, high-value protein
backbone — a single-digit-% ingredient in a general run, or the core of a
bio-specialised mid-training mix. Full TrEMBL is 4–18× larger for mostly
redundant sequence and should be avoided.
