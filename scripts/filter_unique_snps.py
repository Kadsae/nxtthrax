#!/usr/bin/env python3
"""
Remove SNPs from a VCF that are unique to the user's new genome(s)
(i.e. present only in those samples and nowhere else).

Target samples are derived automatically: all .fna/.fasta/.fa files in
the genome directory minus the test genomes listed in exclude_samples.txt.
No manual configuration needed.

Usage:
    python filter_unique_snps.py input.vcf output.vcf
    python filter_unique_snps.py input.vcf output.vcf --genomedir data/genomes --exclude exclude_samples.txt
"""

import argparse
import sys
import vcfpy
from datetime import datetime
from pathlib import Path


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def get_target_samples(genomedir, exclude_file):
    """Return genome filenames in genomedir minus the test genomes in exclude_file."""
    exts = {".fna", ".fasta", ".fa"}
    all_genomes = {p.name for p in Path(genomedir).iterdir() if p.suffix in exts}
    if not all_genomes:
        die(f"No .fna/.fasta/.fa files found in {genomedir}")

    excluded = set()
    if Path(exclude_file).exists():
        with open(exclude_file) as f:
            excluded = {line.strip() for line in f if line.strip()}
    else:
        log(f"WARNING: exclude file not found ({exclude_file}) — treating all genomes as targets")

    targets = all_genomes - excluded
    if not targets:
        die(f"No user genomes left after excluding test genomes from {exclude_file}")
    return targets


def is_mutated(call):
    """Return True if the sample call carries a non-reference allele."""
    if not call or not call.data:
        return False
    gt = call.data.get("GT", "0/0")
    indices = [int(x) for x in gt.replace("|", "/").split("/") if x.isdigit()]
    return sum(indices) > 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf",  help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("--genomedir", "-g", default="data/genomes",
                        help="Directory containing genome files (default: data/genomes)")
    parser.add_argument("--exclude", "-e", default="exclude_samples.txt",
                        help="File listing test genome filenames to exclude (default: exclude_samples.txt)")
    args = parser.parse_args()

    if not Path(args.input_vcf).exists():
        die(f"Input VCF not found: {args.input_vcf}")
    if not Path(args.genomedir).is_dir():
        die(f"Genome directory not found: {args.genomedir}")

    target_samples = get_target_samples(args.genomedir, args.exclude)

    log(f"Input:      {args.input_vcf}")
    log(f"Output:     {args.output_vcf}")
    log(f"Genome dir: {args.genomedir}")
    log(f"Excluded:   {args.exclude}")
    log(f"Targets:    {', '.join(sorted(target_samples))}")

    reader = vcfpy.Reader.from_path(args.input_vcf)
    sample_names = reader.header.samples.names

    # ---------------------------------------------------------------------------
    # sanity check — warn if any genome file has no matching sample in the VCF
    # ---------------------------------------------------------------------------
    missing = target_samples - set(sample_names)
    if missing:
        die(f"Target sample(s) not found in VCF: {', '.join(sorted(missing))}\n"
            "  Make sure genome filenames match the sample names in the VCF header.")

    writer = vcfpy.Writer.from_path(args.output_vcf, reader.header)

    kept = removed = 0

    for record in reader:
        target_mut     = 0
        non_target_mut = 0

        for sample in sample_names:
            call = record.call_for_sample.get(sample)
            if is_mutated(call):
                if sample in target_samples:
                    target_mut += 1
                else:
                    non_target_mut += 1

        # Drop SNPs present only in the user's genomes
        if target_mut > 0 and non_target_mut == 0:
            removed += 1
        else:
            writer.write_record(record)
            kept += 1

    reader.close()
    writer.close()

    log(f"Records kept:    {kept}")
    log(f"Records removed: {removed}")
    log("Done.")


if __name__ == "__main__":
    main()