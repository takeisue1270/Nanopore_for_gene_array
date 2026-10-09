#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""YASS .yop → SVG dot plot (frame-only, height scaled to reference length)

Key options:
  --ref-len N      : force the Y-axis maximum to N (bp). (recommended: REFBASE length)
  --frame-only     : draw only outer frame + alignments (no ticks/labels)

Usage:
  python3 yass_yop_to_svg.py input.yop output.svg --ref-len 1083406 --frame-only
"""

import argparse
import io
import os
import re
import sys
from typing import Tuple

# Alignment line example:
# *(1083406-1083463)(1290780-1290837) Ev: 9.1282 s: 58/58 f
RE_ALN = re.compile(
    r"^\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)\s+Ev:\s+[^s]*\s+s:\s+\d+/\d+\s+([fr])\s*$"
)

def read_bounds(path: str, max_al: int) -> Tuple[int, int, int]:
    """Return (qmax, rmax, n_aln). q=1st coord block, r=2nd coord block."""
    qmax = 0
    rmax = 0
    n = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = RE_ALN.match(line)
            if not m:
                continue
            q1, q2, r1, r2, _ = m.groups()
            q1 = int(q1); q2 = int(q2); r1 = int(r1); r2 = int(r2)
            qmax = max(qmax, q1, q2)
            rmax = max(rmax, r1, r2)
            n += 1
            if n >= max_al:
                break
    return qmax, rmax, n

def svg_header(w: int, h: int) -> str:
    return (f"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\">\n"
            f"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n")

def svg_footer() -> str:
    return "</svg>\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("yop")
    ap.add_argument("svg")
    ap.add_argument("--width", type=int, default=750)
    ap.add_argument("--margin", type=int, default=10)  # small margin for frame-only
    ap.add_argument("--max-al", type=int, default=1_000_000)
    ap.add_argument("--forward-color", default="#000000")
    ap.add_argument("--reverse-color", default="#000000")
    ap.add_argument("--thickness", type=float, default=1.0)
    ap.add_argument("--ref-len", type=int, default=0, help="force Y-axis maximum to this length (bp)")
    ap.add_argument("--frame-only", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.yop):
        sys.stderr.write(f"[ERROR] YOP not found: {args.yop}\n")
        sys.exit(2)

    qmax, rmax, n_aln = read_bounds(args.yop, args.max_al)

    # Empty alignment: still emit a frame
    x_max = qmax if qmax > 0 else 1
    y_max = (args.ref_len if args.ref_len > 0 else rmax) or 1

    fact = x_max / float(args.width) if x_max else 1.0
    if fact <= 0:
        fact = 1.0

    dimX = args.width
    dimY = int(round(y_max / fact))

    totalW = dimX + args.margin * 2
    totalH = dimY + args.margin * 2

    left = args.margin
    top = args.margin
    right = left + dimX
    bottom = top + dimY

    buf = io.StringIO()
    buf.write(svg_header(totalW, totalH))

    # Frame only
    buf.write(f"<rect x=\"{left}\" y=\"{top}\" width=\"{dimX}\" height=\"{dimY}\" fill=\"none\" stroke=\"#000000\"/>\n")

    # Alignments (map: x=query, y=ref; SVG origin is top-left, so invert y)
    fwd = args.forward_color
    rev = args.reverse_color
    thk = args.thickness

    with open(args.yop, "r", encoding="utf-8", errors="ignore") as fh:
        w = buf.write
        for line in fh:
            m = RE_ALN.match(line)
            if not m:
                continue
            q1, q2, r1, r2, fr = m.groups()
            q1 = int(q1); q2 = int(q2); r1 = int(r1); r2 = int(r2)

            x1 = left + int(round(q1 / fact))
            x2 = left + int(round(q2 / fact))

            y1 = top + int(round(dimY - (r1 / fact)))
            y2 = top + int(round(dimY - (r2 / fact)))

            color = fwd if fr == 'f' else rev
            w(f"<line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" stroke=\"{color}\" stroke-width=\"{thk}\"/>\n")

    buf.write(svg_footer())

    with open(args.svg, "w", encoding="utf-8") as out:
        out.write(buf.getvalue())

if __name__ == "__main__":
    main()
