#!/usr/bin/env python3
# Copyright (c) 2026. Released under the Apache License 2.0.
"""
Estimate the token yield of the sequence-first bio corpus under **Marin's
tokenizer** (`marin-community/marin-tokenizer`, which is the Llama-3 128k-vocab
tokenizer), for various database subsets.

It measures per-record token counts on rendered samples (splitting sequence
tokens from metadata tokens), then projects to whole-database totals using
current entry counts.

    # one-time: fetch Marin's tokenizer.json
    hf download marin-community/marin-tokenizer tokenizer.json --local-dir ./marin_tok

    python estimate_tokens.py --tokenizer ./marin_tok/tokenizer.json \
        --sample "Swiss-Prot=sprot_sample.jsonl" \
        --sample "TrEMBL=trembl_sample.jsonl" \
        --sample "Ensembl(pep)=ensembl_pep.sample.jsonl"

Caveats: per-record rates are from samples of ~1k records fetched via the
UniProt search API, which surfaces well-annotated entries first, so measured
TrEMBL *metadata* is an upper bound. The *sequence* term (mean_len x tok/residue)
is organism-independent and robust, so TrEMBL totals are reported as a band.
"""
import argparse
import json
import statistics
from tokenizers import Tokenizer

# Database entry counts. Fetched 2026-07 from the UniProt & UniRef REST APIs
# (x-total-results). Ensembl human is the GRCh38 protein-coding gene count.
DB_COUNTS = {
    "Swiss-Prot (reviewed)": 575_503,
    "TrEMBL (unreviewed)": 149_234_636,
    "UniRef100": 220_919_788,
    "UniRef90": 121_389_642,
    "UniRef50": 38_794_121,
    "Ensembl human (protein-coding genes)": 20_000,
    # new genomic/regulatory record types (human, single species; fetched 2026-07)
    "Ensembl human dogma (coding transcripts)": 20_000,
    "Ensembl human regulatory features": 612_140,
    "Ensembl human splice sites, canonical (donor+acceptor)": 360_000,
}
IDENTITY_LINE_TOKENS = 25  # a minimal ">db:acc name [organism]" + one-sentence lead


def measure(path: str, tok: Tokenizer) -> dict:
    tot, seq_t, slen = [], [], []
    for line in open(path):
        r = json.loads(line)
        tot.append(len(tok.encode(r["text"]).ids))
        seqs = r.get("sequences")
        if seqs:  # multi-sequence (central dogma): count the forms shown in the text
            forms = {k: v for k, v in seqs.items() if k != "cds"}  # cds is a substring of rna; don't double-count
            seq_t.append(sum(len(tok.encode(v).ids) for v in forms.values()))
            slen.append(sum(len(v) for v in forms.values()))
        else:
            seq = r.get("sequence", "")
            seq_t.append(len(tok.encode(seq).ids) if seq else 0)
            slen.append(r.get("seq_len", len(seq)))
    n = len(tot)
    mseq = statistics.mean(seq_t)
    mlen = statistics.mean(slen)
    return {
        "n": n,
        "seq_len": mlen,
        "seq_tok": mseq,
        "meta_tok": statistics.mean(t - s for t, s in zip(tot, seq_t)),
        "total_tok": statistics.mean(tot),
        "median_tok": statistics.median(tot),
        "tok_per_res": (mseq / mlen) if mlen else 0.0,
        "chars_per_tok": (mlen / mseq) if mseq else 0.0,  # for aa records ~ residues/token
    }


def human(n: float) -> str:
    for unit, div in (("T", 1e12), ("B", 1e9), ("M", 1e6), ("K", 1e3)):
        if abs(n) >= div:
            return f"{n / div:.1f}{unit}"
    return f"{n:.0f}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tokenizer", required=True, help="path to tokenizer.json (marin-community/marin-tokenizer)")
    ap.add_argument("--sample", action="append", default=[], metavar="LABEL=PATH", help="repeatable")
    ap.add_argument(
        "--db-mean-aa",
        type=float,
        default=330,
        help="realistic whole-DB mean protein length for the sequence-floor projection "
        "(UniProt ~330-355; the API samples run longer)",
    )
    args = ap.parse_args()

    tok = Tokenizer.from_file(args.tokenizer)
    samples = {}
    for spec in args.sample:
        label, path = spec.split("=", 1)
        samples[label] = measure(path, tok)

    print(f"Tokenizer: {args.tokenizer}  (vocab {tok.get_vocab_size():,})\n")
    print("MEASURED PER-RECORD (sequence-first rendering)")
    hdr = f"{'sample':28} {'n':>5} {'mean aa':>8} {'seq tok':>8} {'meta tok':>9} {'total tok':>10} {'tok/res':>8}"
    print(hdr)
    print("-" * len(hdr))
    for label, s in samples.items():
        print(
            f"{label:28} {s['n']:>5} {s['seq_len']:>8.0f} {s['seq_tok']:>8.0f} "
            f"{s['meta_tok']:>9.0f} {s['total_tok']:>10.0f} {s['tok_per_res']:>8.3f}"
        )

    # ---- projections ----
    def rate_total(label):
        return samples[label]["total_tok"]

    # sequence-token basis: tok/residue is organism-independent (~0.55), so combine
    # it with a realistic whole-DB mean length rather than the (length-skewed) sample.
    seq_basis = "TrEMBL" if "TrEMBL" in samples else next(iter(samples))
    tpr = samples[seq_basis]["tok_per_res"]
    floor_basis = f"{args.db_mean_aa:.0f}aa x {tpr:.2f} + {IDENTITY_LINE_TOKENS}"

    print("\nPROJECTED WHOLE-DATABASE TOKENS (sequence-first)")
    print(f"{'subset':38} {'entries':>13} {'tok/rec':>9} {'total tokens':>14}   basis")
    print("-" * 96)
    rows = []
    if "Swiss-Prot" in samples:
        rows.append(("Swiss-Prot (reviewed)", DB_COUNTS["Swiss-Prot (reviewed)"], rate_total("Swiss-Prot"), "measured (rich)"))
    # TrEMBL as a band: sequence floor -> measured-rich ceiling
    tf = args.db_mean_aa * tpr + IDENTITY_LINE_TOKENS
    rows.append(("TrEMBL (seq + identity line, floor)", DB_COUNTS["TrEMBL (unreviewed)"], tf, floor_basis))
    if "TrEMBL" in samples:
        rows.append(("TrEMBL (rich metadata, ceiling)", DB_COUNTS["TrEMBL (unreviewed)"], rate_total("TrEMBL"), "measured (skewed high)"))
    for ur in ("UniRef100", "UniRef90", "UniRef50"):
        rows.append((f"{ur} (seq + identity line)", DB_COUNTS[ur], tf, floor_basis))
    ens_label = next((k for k in samples if k.startswith("Ensembl")), None)
    if ens_label:
        rows.append(("Ensembl human (protein-coding genes)", DB_COUNTS["Ensembl human (protein-coding genes)"], rate_total(ens_label), f"measured ({ens_label})"))
    # new genomic / regulatory record types
    for label, dbkey in (
        ("Dogma", "Ensembl human dogma (coding transcripts)"),
        ("Regulatory", "Ensembl human regulatory features"),
        ("Splice", "Ensembl human splice sites, canonical (donor+acceptor)"),
    ):
        if label in samples:
            rows.append((dbkey, DB_COUNTS[dbkey], rate_total(label), f"measured ({label})"))

    for name, n, rate, basis in rows:
        print(f"{name:38} {n:>13,} {rate:>9.0f} {human(n * rate):>14}   {basis}")

    print(
        "\nNote: sequence tokens dominate and are robust; TrEMBL/UniRef *metadata* is a band "
        "(the floor assumes a minimal identity line, the ceiling uses annotation-rich samples). "
        "Dedup to UniRef50 is the natural way to cut TrEMBL's redundancy ~4x."
    )


if __name__ == "__main__":
    main()
