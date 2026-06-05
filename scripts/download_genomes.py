#!/usr/bin/env python3
"""
Download B. anthracis genomes from NCBI based on accession numbers in metadata.tsv.

Reads GCF_/GCA_ accessions from the metadata file and downloads the corresponding
genome FASTA files using the NCBI datasets CLI tool.

Usage:
    python download_genomes.py --metadata config/metadata.tsv --outdir data/genomes
    python download_genomes.py --metadata config/metadata.tsv --outdir data/genomes --retry 5
"""

import argparse
import os
import subprocess
import sys
import time
import pandas as pd
from datetime import datetime
from pathlib import Path


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

def check_datasets_tool():
    """Verify the NCBI datasets CLI is available."""
    try:
        r = subprocess.run(["datasets", "--version"], capture_output=True, text=True, timeout=5)
        if r.returncode == 0:
            log(f"NCBI datasets tool found: {r.stdout.strip()}")
            return
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    die("NCBI datasets tool not found.\n"
        "  Install with: conda install -c conda-forge ncbi-datasets-cli\n"
        "  Or see: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download/")


def valid_accession(acc):
    """Return True for GCF_* or GCA_* accessions."""
    if not acc or pd.isna(acc):
        return False
    return str(acc).strip().startswith(("GCF_", "GCA_"))


# ---------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------

def read_metadata(metadata_file):
    """Read TSV and return (dataframe, accession_column_name)."""
    if not os.path.exists(metadata_file):
        die(f"Metadata file not found: {metadata_file}")

    log(f"Reading metadata: {metadata_file}")

    try:
        df = pd.read_csv(metadata_file, sep="\t", dtype=str)
    except Exception as e:
        die(f"Failed to read metadata: {e}")

    log(f"Loaded {len(df)} samples.")

    candidate_cols = [c for c in df.columns
                      if c.lower() in ("accession", "ncbi_accession", "refseq",
                                       "genbank", "assembly_id", "strain")]
    if not candidate_cols:
        die(f"No accession column found. Available columns: {list(df.columns)}\n"
            "  Expected one of: accession, ncbi_accession, refseq, genbank, assembly_id, strain")

    col = candidate_cols[0]
    log(f"Using column '{col}' for accession numbers.")
    return df, col


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------

def download_genome(accession, outdir, retry=3):
    """Download one genome. Returns output path on success, None on failure."""
    accession  = str(accession).strip()
    temp_dir   = os.path.join(outdir, f".temp_{accession}")
    output_fna = os.path.join(outdir, f"{accession}.fna")
    os.makedirs(temp_dir, exist_ok=True)

    for attempt in range(1, retry + 1):
        try:
            # download zip
            r = subprocess.run(
                ["datasets", "download", "genome", "accession", accession,
                 "--filename", os.path.join(temp_dir, f"{accession}.zip"),
                 "--exclude-gff3", "--exclude-protein", "--exclude-seq-report"],
                capture_output=True, text=True, timeout=120
            )
            if r.returncode != 0:
                log(f"  attempt {attempt}/{retry} failed: {r.stderr.strip()}")
                if attempt < retry: time.sleep(5)
                continue

            # extract
            r = subprocess.run(
                ["unzip", "-q", "-o",
                 os.path.join(temp_dir, f"{accession}.zip"), "-d", temp_dir],
                capture_output=True, text=True
            )
            if r.returncode != 0:
                log(f"  attempt {attempt}/{retry} extract failed: {r.stderr.strip()}")
                if attempt < retry: time.sleep(5)
                continue

            # find fasta
            fasta_files = subprocess.run(
                ["find", temp_dir, "-name", "*.fna"],
                capture_output=True, text=True
            ).stdout.strip().splitlines()
            fasta_files = [f for f in fasta_files if f]

            if not fasta_files:
                log(f"  attempt {attempt}/{retry}: no .fna found in archive")
                if attempt < retry: time.sleep(5)
                continue

            subprocess.run(["cp", fasta_files[0], output_fna], check=True)
            subprocess.run(["rm", "-rf", temp_dir], check=True)
            return output_fna

        except subprocess.TimeoutExpired:
            log(f"  attempt {attempt}/{retry} timed out")
            if attempt < retry: time.sleep(5)
        except Exception as e:
            log(f"  attempt {attempt}/{retry} error: {e}")
            if attempt < retry: time.sleep(5)

    subprocess.run(["rm", "-rf", temp_dir], check=False)
    return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", "-m", required=True, help="Path to metadata TSV")
    parser.add_argument("--outdir",   "-o", default="data/genomes",
                        help="Output directory (default: data/genomes)")
    parser.add_argument("--retry",    "-r", type=int, default=3,
                        help="Retry attempts per genome (default: 3)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip genomes already present in output directory")
    args = parser.parse_args()

    check_datasets_tool()
    os.makedirs(args.outdir, exist_ok=True)

    df, acc_col = read_metadata(args.metadata)
    df["accession"] = df[acc_col].astype(str)
    valid = df[df["accession"].apply(valid_accession)].copy()

    if len(valid) == 0:
        die("No valid GCF_/GCA_ accessions found in metadata.")

    skipped_invalid = len(df) - len(valid)
    if skipped_invalid:
        log(f"Skipping {skipped_invalid} sample(s) with missing or invalid accessions.")

    log(f"Downloading {len(valid)} genome(s) to {args.outdir} ...")
    log("--------------------------------------------------------")

    downloaded = failed = skipped = 0

    for i, (_, row) in enumerate(valid.iterrows(), 1):
        acc    = row["accession"]
        strain = row.get("strain", acc)
        outfna = os.path.join(args.outdir, f"{acc}.fna")

        if args.skip_existing and os.path.exists(outfna):
            log(f"[{i}/{len(valid)}] Skipping {acc} (already exists)")
            skipped += 1
            continue

        log(f"[{i}/{len(valid)}] {strain} ({acc})")
        result = download_genome(acc, args.outdir, retry=args.retry)

        if result:
            downloaded += 1
        else:
            log(f"  FAILED: {acc}")
            failed += 1

    # ---------------------------------------------------------------------------
    # summary
    # ---------------------------------------------------------------------------
    log("========================================================"  )
    log("Download complete!")
    log(f"  Downloaded: {downloaded}")
    log(f"  Failed:     {failed}")
    log(f"  Skipped:    {skipped}")
    log("========================================================")

    if failed:
        log("Some downloads failed. Common causes: invalid accession, "
            "NCBI temporarily unavailable, network issues.")
        sys.exit(1)


if __name__ == "__main__":
    main()