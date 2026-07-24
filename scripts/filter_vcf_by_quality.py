#!/usr/bin/env python3
"""
Filter a VCF file by its FILTER column.

Modes
-----
pass-only  (default) Keep only records whose FILTER is exactly PASS.
drop-n               Keep PASS, '.', and any other non-N filter value; drop
                     anything whose FILTER starts with 'N' (N, N:LCB, N:ALN,
                     ... i.e. the parsnp low-confidence flags).
n-only               Keep only records whose FILTER starts with 'N'. This is
                     the exact complement of drop-n, so it is what the default
                     mode throws away -- handy for inspecting the discards.
fail-only            Keep everything that did NOT pass, i.e. the inverse of
                     pass-only.

'.' (no filter applied) is treated as "not assessed", not as "passed". It is
kept in drop-n, and dropped in the other modes unless --keep-missing is given.

Usage:
    python filter_vcf_by_quality.py -i merged.vcf -o pass.vcf
    python filter_vcf_by_quality.py -i merged.vcf -o filtered.vcf --mode drop-n
    python filter_vcf_by_quality.py -i merged.vcf -o discarded.vcf --n-only
    python filter_vcf_by_quality.py -i merged.vcf -o out.vcf --mode fail-only
"""

import argparse
import gzip
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

MODES = ("pass-only", "drop-n", "n-only", "fail-only")

MODE_DESCRIPTIONS = {
    "drop-n": "remove records whose FILTER starts with 'N' (N, N:LCB, N:ALN, ...)",
    "pass-only": "keep only records whose FILTER is exactly PASS",
    "n-only": "keep only records whose FILTER starts with 'N'",
    "fail-only": "keep only records that are not PASS",
}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_out(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def make_predicate(mode, keep_missing=False):
    """Return a function FILTER-value -> bool (True = write the record out)."""

    def normalise(filter_value):
        return filter_value.strip().upper()

    if mode == "drop-n":
        def keep(filter_value):
            v = normalise(filter_value)
            if v in ("PASS", "."):
                return True
            return not v.startswith("N")

    elif mode == "pass-only":
        def keep(filter_value):
            v = normalise(filter_value)
            return v == "PASS" or (keep_missing and v == ".")

    elif mode == "n-only":
        def keep(filter_value):
            v = normalise(filter_value)
            if v == "." and keep_missing:
                return True
            return v.startswith("N")

    elif mode == "fail-only":
        def keep(filter_value):
            v = normalise(filter_value)
            if v == ".":
                return bool(keep_missing)
            return v != "PASS"

    else:
        raise ValueError(f"Unknown mode: {mode}")

    return keep


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", "-i", required=True, help="Input VCF (plain or .gz)")
    parser.add_argument("--output", "-o", required=True, help="Output VCF (.gz writes gzipped)")
    parser.add_argument("--report", "-r", default=None, help="Optional report file")
    parser.add_argument(
        "--keep-missing",
        action="store_true",
        help="Treat FILTER='.' as keepable in pass-only/n-only/fail-only modes",
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--mode", "-m", choices=MODES, default=None,
        help="Filtering mode (default: pass-only)",
    )
    mode_group.add_argument(
        "--pass-only", dest="mode", action="store_const", const="pass-only",
        help="Shorthand for --mode pass-only",
    )
    mode_group.add_argument(
        "--n-only", dest="mode", action="store_const", const="n-only",
        help="Shorthand for --mode n-only",
    )
    mode_group.add_argument(
        "--fail-only", dest="mode", action="store_const", const="fail-only",
        help="Shorthand for --mode fail-only",
    )

    args = parser.parse_args(argv)
    if args.mode is None:
        args.mode = "pass-only"
    return args


def main():
    args = parse_args()

    if not Path(args.input).exists():
        die(f"Input file not found: {args.input}")

    keep = make_predicate(args.mode, keep_missing=args.keep_missing)

    log(f"Input:   {args.input}")
    log(f"Output:  {args.output}")
    log(f"Mode:    {args.mode} -- {MODE_DESCRIPTIONS[args.mode]}")
    if args.keep_missing and args.mode != "drop-n":
        log("         FILTER='.' is being kept (--keep-missing)")

    kept = removed = malformed = 0
    seen_kept = Counter()
    seen_removed = Counter()

    with open_vcf(args.input) as infile, open_out(args.output) as outfile:
        for line in infile:
            # Header lines pass through unchanged
            if line.startswith("#"):
                outfile.write(line)
                continue
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                malformed += 1
                continue

            filter_value = fields[6].strip()
            if keep(filter_value):
                outfile.write(line)
                kept += 1
                seen_kept[filter_value] += 1
            else:
                removed += 1
                seen_removed[filter_value] += 1

    log(f"Records kept:    {kept}")
    log(f"Records removed: {removed}")
    if malformed:
        log(f"Malformed lines skipped (fewer than 7 columns): {malformed}")

    if seen_removed:
        log("Removed by FILTER value:")
        for value, count in seen_removed.most_common():
            log(f"    {value or '<empty>':<20} {count}")

    # ---------------------------------------------------------------------------
    # optional report
    # ---------------------------------------------------------------------------
    if args.report:
        with open(args.report, "w") as f:
            f.write("VCF Filtering Report\n")
            f.write(f"Generated:       {datetime.now().isoformat(timespec='seconds')}\n")
            f.write(f"Input:           {args.input}\n")
            f.write(f"Output:          {args.output}\n")
            f.write(f"Mode:            {args.mode}\n")
            f.write(f"Filter rule:     {MODE_DESCRIPTIONS[args.mode]}\n")
            f.write(f"Keep missing:    {args.keep_missing}\n")
            f.write(f"Records kept:    {kept}\n")
            f.write(f"Records removed: {removed}\n")
            f.write(f"Malformed lines: {malformed}\n")

            if seen_kept:
                f.write("\nKept by FILTER value:\n")
                for value, count in seen_kept.most_common():
                    f.write(f"    {value or '<empty>':<20} {count}\n")
            if seen_removed:
                f.write("\nRemoved by FILTER value:\n")
                for value, count in seen_removed.most_common():
                    f.write(f"    {value or '<empty>':<20} {count}\n")
        log(f"Report written:  {args.report}")

    log("Done.")


if __name__ == "__main__":
    main()
