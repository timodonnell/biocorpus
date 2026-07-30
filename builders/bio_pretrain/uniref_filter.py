#!/usr/bin/env python3
"""Phase-1 of the annotated-UniRef join: subset a UniProtKB flat file to the
records that are UniRef50 cluster *representatives*.

UniRef50 gives the deduplicated protein set (one representative per cluster at
50% identity). Each representative is a real UniProtKB entry, so its annotation
already lives in the flat file (uniprot_{sprot,trembl}.dat.gz). This script does
the cheap half of the join: stream the flat file at byte level, keep only records
whose accession is a UniRef50 representative, and write them round-robin to N
shard files. Each shard is itself a valid UniProt flat file, so phase-2 is just
the ordinary `build_bio_corpus.py uniprot --dat shard.dat ...` (same renderer,
gate, and dedup) run once per shard in parallel.

Representatives are identified from the UniRef50 FASTA headers (RepID); UniParc
(UPI…) representatives are dropped — they are archive-only, with no UniProtKB
annotation.
"""
import argparse
import os
import pickle
import re
import subprocess
import sys

# UniProtKB accession syntax (drops UniParc UPI… representatives).
_UNIPROT_ACC = re.compile(r"^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
_SEP = b"\n//\n"


def _pigz_stdout(path: str):
    return subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 24).stdout


def load_repset(uniref_fasta: str, cache: str | None) -> set:
    if cache and os.path.exists(cache):
        sys.stderr.write(f"[repset] loading cache {cache}\n"); sys.stderr.flush()
        with open(cache, "rb") as f:
            return pickle.load(f)
    sys.stderr.write(f"[repset] scanning {uniref_fasta} headers…\n"); sys.stderr.flush()
    s: set = set()
    out = _pigz_stdout(uniref_fasta)
    for line in out:
        if line[:1] == b">":
            # >UniRef50_Q6GZX4 Cluster… n=… Tax=… TaxID=… RepID=Q6GZX4
            cid = line[1:].split(None, 1)[0].decode()          # UniRef50_Q6GZX4
            acc = cid.split("_", 1)[1] if "_" in cid else cid  # Q6GZX4
            if _UNIPROT_ACC.match(acc):
                s.add(acc)
    if cache:
        with open(cache, "wb") as f:
            pickle.dump(s, f, protocol=5)
        sys.stderr.write(f"[repset] cached {len(s):,} accessions -> {cache}\n"); sys.stderr.flush()
    return s


def stream_records(path: str):
    """Yield raw flat-file record bodies (without the trailing `//`).

    Splits the decompressed stream on the record separator with a single C-level
    ``bytes.split`` per read chunk (keeping only the trailing partial record as
    carry-over). An earlier version sliced the front of the buffer per record
    (``buf = buf[i:]``), which is O(n²) in the chunk size and pegged one core at
    ~2 MB/s; this is linear.
    """
    src = _pigz_stdout(path) if path.endswith(".gz") else open(path, "rb")
    read = src.read
    buf = b""
    for chunk in iter(lambda: read(1 << 23), b""):
        parts = (buf + chunk).split(_SEP)
        buf = parts.pop()          # trailing partial record; completed next round
        yield from parts
    tail = buf.strip()
    if tail:
        yield tail


def accessions(body: bytes) -> list:
    """Accessions on the first AC line (primary + line-1 secondaries).

    A UniRef representative is identified by the entry's *primary* accession, so
    the first AC line is sufficient — and cheap: no regex, no sequence scan.
    Flat-file records always begin with ``ID   ``, so AC is never the first line.
    """
    i = body.find(b"\nAC   ")
    if i < 0:
        return ()
    j = body.find(b"\n", i + 6)
    line = body[i + 6: j] if j >= 0 else body[i + 6:]
    return [tok.strip().decode() for tok in line.split(b";") if tok.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", help="uniprot_{sprot,trembl}.dat[.gz]")
    ap.add_argument("--uniref", required=True, help="uniref50.fasta.gz")
    ap.add_argument("--repset-cache", help="pickle path to cache the representative set")
    ap.add_argument("--out-prefix", help="writes {prefix}.00.dat … {prefix}.NN.dat")
    ap.add_argument("--shards", type=int, default=16)
    ap.add_argument("--repset-only", action="store_true", help="just build the repset cache and exit")
    a = ap.parse_args()

    rep = load_repset(a.uniref, a.repset_cache)
    sys.stderr.write(f"[repset] {len(rep):,} representative accessions\n"); sys.stderr.flush()
    if a.repset_only:
        return
    if not a.inp or not a.out_prefix:
        ap.error("--in and --out-prefix are required unless --repset-only")

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)) or ".", exist_ok=True)
    outs = [open(f"{a.out_prefix}.{i:02d}.dat", "wb", buffering=1 << 20) for i in range(a.shards)]
    seen = kept = k = 0
    for body in stream_records(a.inp):
        seen += 1
        if any(acc in rep for acc in accessions(body)):
            outs[k % a.shards].write(body + _SEP)
            k += 1
            kept += 1
        if seen % 5_000_000 == 0:
            sys.stderr.write(f"  scanned {seen:,}  kept {kept:,}\n"); sys.stderr.flush()
    for o in outs:
        o.close()
    sys.stderr.write(f"[done] scanned {seen:,}  kept {kept:,}  ({100 * kept / max(seen, 1):.1f}%)\n")


if __name__ == "__main__":
    main()
