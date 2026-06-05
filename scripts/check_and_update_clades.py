#!/usr/bin/env python3
"""
check_and_update_clades.py

After bcftools merge, validates that each clade in clades.tsv has a
defining SNP present in the merged VCF.

Logic per clade (site != 0):
  1. If the current site IS in the VCF  -> keep as-is
  2. If the current site is NOT in VCF  -> look up all alternative positions
     for this clade in the Sahl et al. table (NC=7530, Ames coordinates).
     Pick the first position that IS in the VCF, take its ALT allele,
     and update clades.tsv.
  3. If no Sahl alternative found in VCF -> warn and leave unchanged.

Usage (called automatically by run_pipeline.sh):
    python3 scripts/check_and_update_clades.py \
        --vcf  parsnp_new_30.vcf \
        --clades config/clades.tsv \
        --sahl  data/reference/Sahl_et_al_Table3.xlsx
"""

import argparse
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
def get_args():
    p = argparse.ArgumentParser(description="Validate and update clade-defining SNPs.")
    p.add_argument("--vcf",    required=True, help="Merged VCF (parsnp_new_30.vcf)")
    p.add_argument("--clades", required=True, help="config/clades.tsv")
    p.add_argument("--sahl",   required=True, help="Sahl_et_al_Table3.xlsx")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Load VCF into a dict  {position (int): alt_allele (str)}
# ---------------------------------------------------------------------------
def load_vcf(vcf_path):
    vcf_snps = {}
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            pos = int(parts[1])
            alt = parts[4]
            # take only single-nucleotide, non-ambiguous ALTs
            if len(alt) == 1 and alt in "ACGT":
                vcf_snps[pos] = alt
    print(f"  Loaded {len(vcf_snps):,} SNP positions from VCF.")
    return vcf_snps


# ---------------------------------------------------------------------------
# Load Sahl table: returns dict  {clade_name: [list of int positions]}
# Only NC=7530 rows (Ames ancestor coordinate space, same as clades.tsv)
# ---------------------------------------------------------------------------
def load_sahl(sahl_path):
    df = pd.read_excel(sahl_path, sheet_name=0)
    df_ames = df[df["NC"] == 7530][["ID", "Branch"]].dropna()
    df_ames = df_ames[df_ames["ID"].apply(lambda x: str(x).isdigit() or isinstance(x, (int, float)))]
    df_ames["ID"] = df_ames["ID"].astype(int)

    sahl = {}
    for _, row in df_ames.iterrows():
        clade = str(row["Branch"]).strip()
        pos   = int(row["ID"])
        sahl.setdefault(clade, []).append(pos)

    print(f"  Loaded Sahl positions for {len(sahl):,} clade labels (NC=7530).")
    return sahl


# ---------------------------------------------------------------------------
# Normalise clade name for Sahl lookup
# clades.tsv uses names like  A.Br.009_018_(WNA)
# Sahl table uses             A.Br.018 or 009
# We try a few transformations before giving up.
# ---------------------------------------------------------------------------
def sahl_candidates(clade_name):
    """Return a list of candidate strings to look up in the Sahl dict."""
    name = clade_name.strip()
    # strip parenthetical suffix e.g. _(WNA)  _(TEA)
    import re
    base = re.sub(r"_\([^)]*\)$", "", name)   # A.Br.009_018_(WNA) -> A.Br.009_018
    base = re.sub(r"\([^)]*\)$",  "", base).strip()

    candidates = [name, base]

    # handle _018 / /042 combined-clade suffixes
    # e.g. A.Br.009_018 -> try "A.Br.009_018", "A.Br.018 or 009", "A.Br.009"
    m = re.match(r"(A|B|C)\.Br\.(\d+)_(\d+)", base)
    if m:
        prefix, n1, n2 = m.group(1), m.group(2), m.group(3)
        candidates += [
            f"{prefix}.Br.{n1} or {n2}",
            f"{prefix}.Br.{n2} or {n1}",
            f"{prefix}.Br.{n1}",
            f"{prefix}.Br.{n2}",
        ]

    m2 = re.match(r"(A|B|C)\.Br\.(\d+)/(\d+)", base)
    if m2:
        prefix, n1, n2 = m2.group(1), m2.group(2), m2.group(3)
        candidates += [
            f"{prefix}.Br.{n1}/{n2}",
            f"{prefix}.Br.{n1}",
            f"{prefix}.Br.{n2}",
        ]

    # deduplicate while preserving order
    seen, out = set(), []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = get_args()

    print("\n[check_and_update_clades] Loading inputs...")
    vcf_snps = load_vcf(args.vcf)
    sahl     = load_sahl(args.sahl)
    clades   = pd.read_csv(args.clades, sep="\t")

    # clades.tsv columns: clade  gene  site  alt
    # site==0  means "placeholder / not defined" -> skip those
    updated_rows  = []
    missing_rows  = []
    ok_rows       = []

    print("\n[check_and_update_clades] Checking each clade...\n")

    for _, row in clades.iterrows():
        clade = str(row["clade"]).strip()
        site  = int(row["site"])
        alt   = str(row["alt"]).strip()

        # placeholders - nothing to validate
        if site == 0:
            ok_rows.append(clade)
            continue

        # --- Case 1: current SNP is in the VCF ---
        if site in vcf_snps:
            ok_rows.append(clade)
            continue

        # --- Case 2: current SNP missing - search Sahl alternatives ---
        candidates = sahl_candidates(clade)
        sahl_positions = []
        matched_key = None
        for key in candidates:
            if key in sahl:
                sahl_positions = sahl[key]
                matched_key = key
                break

        if not sahl_positions:
            missing_rows.append({
                "clade": clade, "old_site": site,
                "reason": "clade not found in Sahl table"
            })
            print(f"  WARN  {clade:40s}  site {site:>8d}  not in VCF. Clade absent from Sahl table")
            continue

        # find first Sahl position that IS in the VCF
        replacement = None
        for sahl_pos in sahl_positions:
            if sahl_pos in vcf_snps:
                replacement = sahl_pos
                break

        if replacement is None:
            missing_rows.append({
                "clade": clade, "old_site": site,
                "reason": f"no Sahl position present in VCF (checked {len(sahl_positions)} positions)"
            })
            print(f"  WARN  {clade:40s}  site {site:>8d}  not in VCF! No Sahl alternative found in VCF either")
            continue

        new_alt = vcf_snps[replacement]
        updated_rows.append({
            "clade": clade,
            "old_site": site, "old_alt": alt,
            "new_site": replacement, "new_alt": new_alt,
            "sahl_key": matched_key
        })
        print(f"  UPDATE {clade:40s}  {site:>8d}/{alt} -> {replacement:>8d}/{new_alt}  (Sahl key: {matched_key})")

        # patch the row in the dataframe
        idx = clades.index[clades["clade"] == row["clade"]].tolist()
        if idx:
            clades.at[idx[0], "site"] = replacement
            clades.at[idx[0], "alt"]  = new_alt

    # --- Summary ---
    print(f"\n[check_and_update_clades] Summary")
    print(f"  OK (present in VCF or placeholder): {len(ok_rows)}")
    print(f"  Updated with Sahl alternative:      {len(updated_rows)}")
    print(f"  Could not resolve (warnings above): {len(missing_rows)}")

    # --- Write updated clades.tsv (only if something changed) ---
    if updated_rows:
        clades.to_csv(args.clades, sep="\t", index=False)
        print(f"\n  Clades file updated: {args.clades}")
    else:
        print(f"\n  No changes needed! Clades file unchanged.")

    # --- Write a log of changes ---
    log_path = args.clades.replace(".tsv", "_update_log.tsv")
    if updated_rows or missing_rows:
        log_df = pd.DataFrame(
            updated_rows + [
                {"clade": r["clade"], "old_site": r["old_site"], "old_alt": "",
                 "new_site": "", "new_alt": "", "sahl_key": r["reason"]}
                for r in missing_rows
            ]
        )
        log_df.to_csv(log_path, sep="\t", index=False)
        print(f"  Change log written: {log_path}")

    print()

    # exit non-zero only if there are unresolvable clades so the pipeline can decide
    if missing_rows:
        sys.exit(2)


if __name__ == "__main__":
    main()
