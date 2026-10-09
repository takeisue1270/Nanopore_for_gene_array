#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract reads that map to TARGET (tseqname) from merged.*.fastq using mappy.
Output: <input>.tgt.fasta
ID: keep full ID from '@' to before first TAB (if present).
"""

import argparse
import sys
from typing import Iterator, Tuple, Optional

try:
    import mappy as mp
except ImportError:
    sys.stderr.write("[ERROR] mappy not found. Please `pip install mappy`.\n")
    sys.exit(1)


def parse_fastq(path: str) -> Iterator[Tuple[str, str]]:
    with open(path, "r", newline="") as fh:
        while True:
            id_line = fh.readline()
            if not id_line:
                break
            seq_line = fh.readline()
            plus_line = fh.readline()
            qual_line = fh.readline()
            if not (seq_line and plus_line and qual_line):
                break
            if not id_line.startswith("@"):
                continue

            # keep full ID from '@' to before first TAB
            rid = id_line.rstrip("\n").split("\t", 1)[0]
            seq = seq_line.strip()
            yield rid, seq


def has_hit(aligner: mp.Aligner, seq: str, min_mapq: int) -> bool:
    best_mapq = -1
    for a in aligner.map(seq):
        if a.q_en > a.q_st and a.mapq > best_mapq:
            best_mapq = a.mapq
    return best_mapq >= min_mapq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tseqname", help="target fasta")
    ap.add_argument("fastqname", help="input fastq (merged.*.fastq)")
    ap.add_argument("--threads", type=int, default=14)
    ap.add_argument("--preset", default="map-ont")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--out", default=None, help="output fasta (default: <fastq>.tgt.fasta)")
    args = ap.parse_args()

    out_fa = args.out or (args.fastqname.rsplit(".fastq", 1)[0] + ".tgt.fasta")

    tgt_aln = mp.Aligner(args.tseqname, preset=args.preset, n_threads=args.threads)
    if not tgt_aln:
        sys.stderr.write(f"[ERROR] failed to init aligner: {args.tseqname}\n")
        sys.exit(1)

    n_in = 0
    n_out = 0
    with open(out_fa, "w") as out:
        for rid, seq in parse_fastq(args.fastqname):
            n_in += 1
            if has_hit(tgt_aln, seq, args.min_mapq):
                out.write(f">{rid}\n{seq}\n")
                n_out += 1

    sys.stderr.write(f"[INFO] processed={n_in}, passed(target)={n_out}, output={out_fa}\n")


if __name__ == "__main__":
    main()
