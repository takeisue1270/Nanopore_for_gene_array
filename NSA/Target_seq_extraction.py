# -*- coding: utf-8 -*-
# Fast version of Target_seq_extraction.py (mappy / minimap2 Python API)
#
# Usage:
#   python3 Target_seq_extraction.py flanking_up.fa flanking_down.fa target.fa input.fastq
# Options:
#   --threads 14 --min-mapq 20 --preset map-ont

import sys
import re
import argparse
from typing import Iterator, Tuple, Optional

try:
    import mappy as mp
except ImportError:
    sys.stderr.write(
        "[ERROR] mappy was not found. Please run `pip install mappy` (Python bindings for minimap2).\n"
    )
    sys.exit(1)


def parse_fastq(path: str) -> Iterator[Tuple[str, str]]:
    """
    Simple FASTQ parser that reads records in 4-line blocks and returns
    (read_id, sequence).
    """
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
                # Skip unexpected formats
                continue
            rid = id_line.strip().split()[0]  # use only the first field
            seq = seq_line.strip()
            yield (rid, seq)


def best_hit_q_coords(
    aligner: mp.Aligner, seq: str, min_mapq: int
) -> Optional[Tuple[int, int]]:
    """
    Align the given read sequence to the provided aligner and return
    the aligned query interval [q_st, q_en) with the highest MAPQ
    passing the given threshold. Returns None if no valid hit is found.
    """
    best = None
    best_mapq = -1
    for a in aligner.map(seq):
        # a: AlignedSegment-like object with attributes:
        # ctg, r_st, r_en, strand, mapq, q_st, q_en, cigar, ...
        if a.mapq >= min_mapq and a.q_en > a.q_st:
            if a.mapq > best_mapq:
                best_mapq = a.mapq
                best = (a.q_st, a.q_en)
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument("fseq1name")   # upstream flanking FASTA
    p.add_argument("fseq2name")   # downstream flanking FASTA
    p.add_argument("tseqname")    # target FASTA
    p.add_argument("fastqname")   # input FASTQ
    p.add_argument("--threads", type=int, default=14)
    p.add_argument("--min-mapq", type=int, default=10)
    p.add_argument(
        "--preset", default="map-ont"
    )  # ONT preset (modify if needed)
    args = p.parse_args()

    # Output file name (same convention as the original script)
    outputname = re.sub(r"\.fastq$", r".target.positive", args.fastqname)

    # Initialize aligners (mappy builds indices from reference FASTA files)
    # References are small, so loading overhead is minimal
    up_aln = mp.Aligner(args.fseq1name, preset=args.preset, n_threads=args.threads)
    if not up_aln:
        sys.stderr.write(f"[ERROR] Failed to initialize aligner: {args.fseq1name}\n")
        sys.exit(1)

    down_aln = mp.Aligner(args.fseq2name, preset=args.preset, n_threads=args.threads)
    if not down_aln:
        sys.stderr.write(f"[ERROR] Failed to initialize aligner: {args.fseq2name}\n")
        sys.exit(1)

    tgt_aln = mp.Aligner(args.tseqname, preset=args.preset, n_threads=args.threads)
    if not tgt_aln:
        sys.stderr.write(f"[ERROR] Failed to initialize aligner: {args.tseqname}\n")
        sys.exit(1)

    # Write header once
    with open(outputname, "w") as out:
        out.write("Repeatsize\tReadID\tSequence\n")

    # Process each read
    count_in = 0
    count_pass = 0
    with open(outputname, "a") as out:
        for rid, seq in parse_fastq(args.fastqname):
            count_in += 1

            up = best_hit_q_coords(up_aln, seq, args.min_mapq)
            if not up:
                continue
            down = best_hit_q_coords(down_aln, seq, args.min_mapq)
            if not down:
                continue
            tgt = best_hit_q_coords(tgt_aln, seq, args.min_mapq)
            if not tgt:
                continue

            # Estimate repeat size from read coordinates where up/down flanks align
            u0, u1 = up
            d0, d1 = down
            candidates = (
                abs(u0 - d0),
                abs(u0 - d1),
                abs(u1 - d0),
                abs(u1 - d1),
            )
            repeatsize = min(candidates)

            out.write(f"{repeatsize}\t{rid}\t{seq}\n")
            count_pass += 1

    sys.stderr.write(
        f"[INFO] processed={count_in}, passed(up/down/target)={count_pass}, "
        f"output={outputname}\n"
    )


if __name__ == "__main__":
    main()
