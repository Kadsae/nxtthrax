#!/usr/bin/env python3
"""
Filter a VCF file by removing records whose FILTER value starts with 'N'.

This covers N, N:LCB, N:ALN, and any other N-prefixed parsnp filter flags,
which all indicate low-quality or low-confidence variant calls.

Usage:
    python filter_vcf_by_quality.py --input merged.vcf --output filtered.vcf
"""

import argparse
import gzip
import sys
from datetime import datetime
from pathlib import Path


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def should_keep(filter_value):
    """Keep PASS and '.' unconditionally; drop anything that starts with 'N'."""
    v = filter_value.strip().upper()
    if v in ("PASS", "."):
        return True
    return not v.startswith("N")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  "-i", required=True, help="Input VCF (plain or .gz)")
    parser.add_argument("--output", "-o", required=True, help="Output VCF")
    parser.add_argument("--report", "-r", default=None,  help="Optional report file")
    args = parser.parse_args()

    if not Path(args.input).exists():
        die(f"Input file not found: {args.input}")

    log(f"Input:   {args.input}")
    log(f"Output:  {args.output}")
    log("Removing records with FILTER starting with 'N' (N, N:LCB, N:ALN, ...)")

    kept = removed = 0

    with open_vcf(args.input) as infile, open(args.output, "w") as outfile:
        for line in infile:
            # Header lines pass through unchanged
            if line.startswith("#"):
                outfile.write(line)
                continue
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                removed += 1
                continue

            if should_keep(fields[6]):
                outfile.write(line)
                kept += 1
            else:
                removed += 1

    log(f"Records kept:    {kept}")
    log(f"Records removed: {removed}")

    # ---------------------------------------------------------------------------
    # optional report
    # ---------------------------------------------------------------------------
    if args.report:
        with open(args.report, "w") as f:
            f.write("VCF Filtering Report\n")
            f.write(f"Input:           {args.input}\n")
            f.write(f"Output:          {args.output}\n")
            f.write("Filter rule:     remove FILTER starting with 'N'\n")
            f.write(f"Records kept:    {kept}\n")
            f.write(f"Records removed: {removed}\n")
        log(f"Report written:  {args.report}")

    log("Done.")


if __name__ == "__main__":
    main()