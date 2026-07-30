# biocorpus

Explorations in biology-focused pre-training datasets for LLMs.

## biocorpus — a sequence-first biology pretraining corpus

**21.1 M documents · ~9.15 B tokens** (Llama-3 / `marin-tokenizer`), built from UniProt and
Ensembl. Every document places a biological **sequence first**, then its real annotation —
training a model to read structure/function *from* sequence. Deduplicated at UniRef50,
filtered to genuinely learnable annotation (no uncharacterized filler, no markup tags), with
**verified** DNA→RNA→protein records.

- 📦 **Dataset:** [huggingface.co/buckets/timodonnell/biocorpus](https://huggingface.co/buckets/timodonnell/biocorpus)
- 📊 **Report:** [timodonnell.github.io/biocorpus/analyses/biocorpus-v2-report.html](https://timodonnell.github.io/biocorpus/analyses/biocorpus-v2-report.html)
- 🛠️ **Builder:** [`builders/bio_pretrain/`](builders/bio_pretrain/) — one `biopython`-only script per source, plus the annotated-UniRef join

| source | documents | tokens |
|---|--:|--:|
| TrEMBL representatives (annotated, UniRef50-deduped) | 20,091,178 | 8.49 B |
| Swiss-Prot (both orderings) | 950,174 | 0.58 B |
| human central-dogma + splice | 96,768 | 0.07 B |

## Analyses

- [Review of trillionlabs/TheBioCollection](https://timodonnell.github.io/biocorpus/analyses/thebiocollection-review.html) — the shortcomings this corpus was built to fix (no provenance, empty-slot template bugs, unlearnable "verifiable" tasks, DNA/RNA/protein never co-occurring).
