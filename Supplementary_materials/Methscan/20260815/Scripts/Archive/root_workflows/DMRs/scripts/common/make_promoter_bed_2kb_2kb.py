#!/usr/bin/env python3
"""Generate a strand-aware GENCODE promoter BED (TSS ±2 kb)."""

import argparse
import re
from pathlib import Path


def parse_attributes(text):
    attrs = {}
    for item in text.strip().split(";"):
        match = re.match(r'(\S+)\s+"([^"]+)"', item.strip())
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--upstream", type=int, default=2000)
    parser.add_argument("--downstream", type=int, default=2000)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with Path(args.gtf).open() as source, output.open("w") as target:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attr = fields
            attrs = parse_attributes(attr)
            start, end = int(start), int(end)
            tss = start if strand == "+" else end if strand == "-" else None
            if tss is None:
                continue
            # GTF: 1-based inclusive; BED: 0-based half-open.
            bed_start = max(0, tss - args.upstream - 1)
            bed_end = tss + args.downstream
            target.write("\t".join(map(str, [
                chrom, bed_start, bed_end,
                attrs.get("gene_id", "NA").split(".")[0],
                attrs.get("gene_name", attrs.get("gene_id", "NA")),
                attrs.get("gene_type", attrs.get("gene_biotype", "NA")),
                strand, tss,
            ])) + "\n")
            n += 1
    print(f"Written: {output}")
    print(f"Promoter definition: TSS - {args.upstream} bp to TSS + {args.downstream} bp")
    print(f"Number of promoter regions: {n}")


if __name__ == "__main__":
    main()
