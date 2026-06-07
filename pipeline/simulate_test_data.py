#!/usr/bin/env python
"""
simulate_test_data.py — tiny synthetic Nanopore amplicon dataset for end-to-end testing.

Builds a MinKNOW-style fastq_pass/barcodeNN/ layout with amplicon reads (including
ORIGIN-SPANNING amp3 reads) for a few HPV types, plus one co-infection barcode. Reads
carry a realistic `basecall_model_version_id=...` header so model auto-selection can be
tested, and a small substitution error rate so consensus is still recoverable.

Reads are drawn from the CANONICAL bundle references; the amp3 read is the wrapped insert
(canonical[L_start:] + canonical[:R_end]) so it genuinely spans the circular origin.
"""
import argparse, gzip, os, random, re
from pathlib import Path
import hpvpipe as H

# 400bps (matches the installed r1041_e82_400bps_{fast,hac,sup} models — real GridION /
# PromethION P2 Solo output is 400bps; the old 260bps id selected a non-existent model).
MODEL_HDR = "basecall_model_version_id=r1041_e82_400bps_hac_v520"

def mutate(seq, rate, rng):
    out = []
    for b in seq:
        if rng.random() < rate:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        else: out.append(b)
    return "".join(out)

def amplicon_span_seq(canon, a, L):
    """full amplicon span LEFT.start..RIGHT.end (wrapped for the junction amplicon)."""
    s, e = a["L"][0], a["R"][1]
    return canon[s:e] if s <= e else canon[s:] + canon[:e]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme-dir", required=True)
    ap.add_argument("--out", default="fastq_pass")
    ap.add_argument("--depth", type=int, default=40)
    ap.add_argument("--error", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    # barcode -> list of types (barcode03 = co-infection)
    layout = {"barcode01": ["HPV16"], "barcode02": ["HPV58"], "barcode03": ["HPV16", "HPV58"]}
    samples = {"barcode01": "SAMPLE_HPV16", "barcode02": "SAMPLE_HPV58", "barcode03": "SAMPLE_COINF"}
    outroot = Path(a.out); outroot.mkdir(parents=True, exist_ok=True)

    for bc, types in layout.items():
        d = outroot / bc; d.mkdir(parents=True, exist_ok=True)
        reads = []
        for t in types:
            ref_fa, primer_bed, _ = H.scheme_files(a.scheme_dir, t)
            canon = next(iter(H.read_fasta(ref_fa).values())); L = len(canon)
            amps = H.load_amplicons(primer_bed); jamp = H.junction_amp(amps, L)
            for amp, av in amps.items():
                if "L" not in av or "R" not in av: continue
                span = amplicon_span_seq(canon, av, L)
                tag = "JUNCTION" if amp == jamp else "amp"
                for i in range(a.depth):
                    s = mutate(span, a.error, rng)
                    if rng.random() < 0.5: s = H.revcomp(s)     # random strand
                    rid = f"{t}_{amp}_{tag}_{i} {MODEL_HDR}"
                    reads.append((rid, s))
        rng.shuffle(reads)
        with gzip.open(d / f"{bc}_reads_0.fastq.gz", "wt") as fh:
            for rid, s in reads:
                fh.write(f"@{rid}\n{s}\n+\n{'I'*len(s)}\n")
        print(f"  {bc} ({'+'.join(types)}): {len(reads)} reads -> {d}")

    # sample sheet
    ss = outroot.parent / "sample_sheet.csv"
    with open(ss, "w") as fh:
        fh.write("barcode,sample_id,expected_types,notes\n")
        for bc, types in layout.items():
            fh.write(f"{bc},{samples[bc]},{';'.join(types)},synthetic\n")
    print(f"  sample sheet -> {ss}")

if __name__ == "__main__":
    main()
