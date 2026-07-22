#!/usr/bin/env python3
# Copyright (c) 2026. Released under the Apache License 2.0.
"""
Cross-file exact-sequence dedup for corpus JSONL shards.

Keeps the **first** occurrence of each sequence; pass files in **priority order**
(richest first, e.g. Swiss-Prot before UniRef before Ensembl). Keys are
namespaced so like dedups against like:

  * protein records (`seq_type=aa`)  -> dedup on the amino-acid sequence
  * DNA records (`seq_type=dna`)      -> dedup on the DNA window
  * RNA records (`seq_type=rna`)      -> dedup on the RNA
  * central-dogma records            -> dedup by `id` only (never dropped
                                         against a protein/DNA file — the value
                                         is the DNA->RNA->protein mapping)
  * records with no sequence          -> always kept

So the same human protein coming from Swiss-Prot, its UniRef50 representative,
and an Ensembl peptide collapses to one; a dogma record whose protein also
appears in Swiss-Prot is untouched.

    python dedup.py sprot.jsonl uniref50.jsonl ensembl_pep.jsonl dogma_*.jsonl --out corpus.jsonl

Memory scales with the number of *unique* records (a 16-byte hash each: ~1.6 GB
per 100 M). Inputs/outputs may be plain or `.gz`.
"""
import argparse
import collections
import gzip
import hashlib
import io
import json
import sys


def _open(path: str, mode: str = "rt") -> io.TextIOBase:
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _key(rec: dict):
    """16-byte namespaced dedup key, or None to always keep."""
    if rec.get("entity_type") == "central_dogma":
        payload = b"D\x00" + (rec.get("id") or "").encode("utf-8")
    else:
        seq = rec.get("sequence") or ""
        if not seq:
            return None
        payload = (rec.get("seq_type") or "?").encode("ascii") + b"\x00" + seq.encode("utf-8")
    return hashlib.blake2b(payload, digest_size=16).digest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="JSONL[.gz] shards, in priority order (first wins)")
    ap.add_argument("--out", required=True, help="deduped output JSONL[.gz]")
    args = ap.parse_args()

    seen: set = set()
    kept = dropped = 0
    kept_by = collections.Counter()
    dropped_by = collections.Counter()

    with _open(args.out, "wt") as out:
        for path in args.inputs:
            k0, d0 = kept, dropped
            for line in _open(path):
                if not line.strip():
                    continue
                rec = json.loads(line)
                src = rec.get("source", "?")
                key = _key(rec)
                if key is not None and key in seen:
                    dropped += 1
                    dropped_by[src] += 1
                    continue
                if key is not None:
                    seen.add(key)
                out.write(line if line.endswith("\n") else line + "\n")
                kept += 1
                kept_by[src] += 1
            print(f"  {path}: kept {kept - k0:,}  dropped {dropped - d0:,}", file=sys.stderr)

    print(f"\ntotal kept {kept:,}  dropped {dropped:,}  ({100 * dropped / max(kept + dropped, 1):.1f}% redundant)", file=sys.stderr)
    for src in sorted(set(kept_by) | set(dropped_by)):
        print(f"    {src:22} kept {kept_by[src]:>12,}  dropped {dropped_by[src]:>12,}", file=sys.stderr)


if __name__ == "__main__":
    main()
