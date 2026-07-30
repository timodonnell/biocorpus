#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# biocorpus v2 — full build under the annotation-quality rules
#   * no <tag> markup; only records carrying real, learnable annotation (gate)
#   * protein bulk is the DEDUPLICATED + ANNOTATED UniRef50 set, obtained by a
#     local flat-file join (no REST): keep the UniProtKB entry of each UniRef50
#     representative, render its annotation. Swiss-Prot is added in full; human
#     central-dogma / splice records give the genomic accent.
#
# Inputs (already downloaded into $DL):
#   sprot.dat.gz  trembl.dat.gz  uniref50.fasta.gz  repset50.pkl
#   human.{dna.fa,gff3.gz,cdna.fa.gz,cds.fa.gz,pep.fa.gz}
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

DL=/data/biocorpus_out/dl
OUT=/data/biocorpus_out/v2
SHARDS=32
mkdir -p "$OUT"
run() { echo "[$(date +%H:%M:%S)] $*"; }

# --- [A] Swiss-Prot: all annotated entries, both orderings (recognition+design) ---
if [ ! -s "$OUT/sprot.jsonl" ]; then
  run "[A] Swiss-Prot (both orderings)"
  python3 build_bio_corpus.py uniprot --dat "$DL/sprot.dat.gz" \
    --ordering both --limit 0 --out "$OUT/sprot.jsonl" > "$OUT/log.sprot" 2>&1 &
fi
SPROT_PID=${!:-}

# --- [B1] TrEMBL -> UniRef50-representative flat-file shards (phase-1 filter) ---
if ! ls "$OUT"/trembl_rep.*.dat >/dev/null 2>&1; then
  run "[B1] filtering TrEMBL to UniRef50 representatives -> $SHARDS shards"
  python3 uniref_filter.py --in "$DL/trembl.dat.gz" --uniref "$DL/uniref50.fasta.gz" \
    --repset-cache "$DL/repset50.pkl" --out-prefix "$OUT/trembl_rep" --shards "$SHARDS" \
    2>&1 | tee "$OUT/log.filter"
fi

# --- [B2] build each rep shard in parallel (phase-2: existing renderer + gate) ---
run "[B2] building $SHARDS rep shards (sequence_first)"
pids=()
for f in "$OUT"/trembl_rep.*.dat; do
  b=$(basename "$f" .dat)
  [ -s "$OUT/$b.jsonl" ] && continue
  python3 build_bio_corpus.py uniprot --dat "$f" \
    --ordering sequence_first --limit 0 --out "$OUT/$b.jsonl" > "$OUT/log.$b" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]:-}"; do wait "$p" || true; done

# --- [C] Human central-dogma + splice (local genome, offline) ---
if [ ! -s "$OUT/dogma_human.jsonl" ]; then
  run "[C] human central-dogma"
  python3 build_bio_corpus.py ensembl-dogma --gff "$DL/human.gff3.gz" \
    --genome "$DL/human.dna.fa" --cdna "$DL/human.cdna.fa.gz" \
    --cds "$DL/human.cds.fa.gz" --pep "$DL/human.pep.fa.gz" \
    --view all --min-exons 2 --limit 0 --out "$OUT/dogma_human.jsonl" > "$OUT/log.dogma" 2>&1 || true
fi
if [ ! -s "$OUT/splice_human.jsonl" ]; then
  run "[C] human splice junctions"
  python3 build_bio_corpus.py ensembl-splice --gff "$DL/human.gff3.gz" \
    --genome "$DL/human.dna.fa" --site both --limit 0 \
    --out "$OUT/splice_human.jsonl" > "$OUT/log.splice" 2>&1 || true
fi

[ -n "$SPROT_PID" ] && wait "$SPROT_PID" || true

# --- [D] cross-file dedup (priority: Swiss-Prot > TrEMBL reps > genomic) ---
run "[D] cross-file dedup"
python3 dedup.py \
  "$OUT/sprot.jsonl" "$OUT"/trembl_rep.*.jsonl \
  "$OUT/dogma_human.jsonl" "$OUT/splice_human.jsonl" \
  --out "$OUT/corpus_v2.jsonl" 2>&1 | tee "$OUT/log.dedup"

run "done -> $OUT/corpus_v2.jsonl  ($(wc -l < "$OUT/corpus_v2.jsonl") docs)"
