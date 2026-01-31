#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract upstream/downstream flanks around the best target (tseqname) hit using mappy.

Goal:
- From merged.*.fastq, find reads that align to target.fa (tseqname)
- For each read with a target hit (MAPQ filtered), extract:
    upstream = read[0:q_st]
    downstream = read[q_en:]
- Write them to a FASTA file next to the input fastq.

Usage:
  python3 Extract_tgt_flanks_mappy.py target.fa input.fastq \
      --threads 14 --min-mapq 10 --min-flank-len 1

Output (default):
  <input>.tgt_flanks.fa
  headers: >readID_up, >readID_down
"""

import argparse
import re
import sys
from typing import Iterator, Tuple, Optional

try:
    import mappy as mp
except ImportError:
    sys.stderr.write("[ERROR] mappy not found. Install with: pip install mappy\n")
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
            rid = id_line.strip().split()[0]
            seq = seq_line.strip()
            yield rid, seq

def best_hit_q_coords(aligner: mp.Aligner, seq: str, min_mapq: int, min_span: int) -> Optional[Tuple[int, int]]:
    best = None
    best_score = (-1, -1)  # (mapq, span)
    for a in aligner.map(seq):
        span = a.q_en - a.q_st
        if a.mapq >= min_mapq and span >= min_span:
            score = (a.mapq, span)
            if score > best_score:
                best_score = score
                best = (a.q_st, a.q_en)
    return best

def main():
    p = argparse.ArgumentParser()
    p.add_argument("tseqname", help="target fasta (tseqname)")
    p.add_argument("fastqname", help="input FASTQ (e.g., merged.*.fastq)")
    p.add_argument("--threads", type=int, default=14)
    p.add_argument("--min-mapq", type=int, default=10)
    p.add_argument("--min-span", type=int, default=50, help="minimum aligned span on read")
    p.add_argument("--min-flank-len", type=int, default=1, help="minimum flank length to output")
    p.add_argument("--preset", default="map-ont")
    p.add_argument("--out", default="", help="output fasta path (default: <input>.tgt_flanks.fa)")
    args = p.parse_args()

    out = args.out or re.sub(r"\.fastq$", ".tgt_flanks.fa", args.fastqname)
    if out == args.fastqname:
        out = args.fastqname + ".tgt_flanks.fa"

    aln = mp.Aligner(args.tseqname, preset=args.preset, n_threads=args.threads, best_n=5)
    if not aln:
        sys.stderr.write(f"[ERROR] failed to init aligner: {args.tseqname}\n")
        sys.exit(2)

    n_in = n_hit = n_out = 0
    with open(out, "w") as fo:
        for rid, seq in parse_fastq(args.fastqname):
            n_in += 1
            hit = best_hit_q_coords(aln, seq, args.min_mapq, args.min_span)
            if not hit:
                continue
            n_hit += 1
            q_st, q_en = hit
            up = seq[:q_st]
            dn = seq[q_en:]
            if len(up) >= args.min_flank_len:
                fo.write(f">{rid}_up\n{up}\n")
                n_out += 1
            if len(dn) >= args.min_flank_len:
                fo.write(f">{rid}_down\n{dn}\n")
                n_out += 1

    sys.stderr.write(f"[INFO] processed={n_in}, tgt_hit={n_hit}, written={n_out}, output={out}\n")

if __name__ == "__main__":
    main()
