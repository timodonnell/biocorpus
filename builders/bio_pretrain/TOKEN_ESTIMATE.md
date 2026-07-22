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

| sample | n | mean len | seq tok | meta tok | total tok |
|---|--:|--:|--:|--:|--:|
| Swiss-Prot (reviewed) | 1000 | 555 aa | 306 | 732 | 1038 |
| TrEMBL (unreviewed) | 1000 | 779 aa | 429 | 586 | 1014 |
| Ensembl human (pep) | 30 | 351 aa | 194 | 140 | 333 |
| **Ensembl dogma** (DNA+RNA+protein) | 20 | 3209 | 1664 | 141 | 1805 |
| **Ensembl regulatory** feature | 25 | 474 nt | 242 | 90 | 331 |
| **Ensembl splice** site (donor/acceptor) | 30 | 200 nt | 104 | 134 | 239 |

(Dogma "mean len" sums the DNA + RNA + protein forms; its `seq tok` excludes the
CDS, which is a substring of the RNA.)

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
| **Ensembl dogma** (DNA+RNA+protein, human coding tx) | ~20,000 | 1,805 | **~36 M** | measured |
| **Ensembl regulatory** features (human) | 612,140 | 331 | **~203 M** | measured |
| **Ensembl splice** sites (human canonical, donor+acceptor) | ~360,000 | 239 | **~86 M** | measured |

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
- **The new genomic/regulatory record types are cheap** (human, single species):
  dogma ~36 M, regulatory features ~203 M, splice sites ~86 M — together
  **~0.3 B tokens**, a rounding error next to the protein backbone, and mostly
  *novel* signal (transcription/splicing/translation, cis-regulatory elements)
  absent from UniProt. They're DNA-heavy (~0.52 tok/nt). Going multi-species
  (Ensembl has ~300 vertebrates; Ensembl Genomes thousands more) scales these
  ~100–300×, at which point dedup/selection matters — but the **human set alone
  (~0.3 B) is an easy, high-value include**.
- **Sequences are expensive tokens.** ~50–90% of every record here is sequence
  tokens carrying little signal for a *text* model. That is the strongest argument
  for (a) keeping the bio fraction small, (b) deduping hard (UniRef50), and
  (c) considering sequence-aware tokenization if the fraction grows.

**Rule of thumb:** Swiss-Prot (~0.6 B, both orderings ≈ 1.2 B) + UniRef50 (~8 B)
+ curated Ensembl gene/annotation records is a ~10 B-token, high-value protein
backbone — a single-digit-% ingredient in a general run, or the core of a
bio-specialised mid-training mix. Full TrEMBL is 4–18× larger for mostly
redundant sequence and should be avoided.
