#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG - only change these if you re-organize the folder layout
# ---------------------------------------------------------------------------
GENOME_DIR="data/genomes"                           # all genomes to include in the analysis
REFERENCE="data/reference/Ames-Ancestor.fna"       # parsnp reference genome
PARSNP_OUTDIR="01_parsnp_out_new"                  # fresh parsnp output folder
FILTERED_VCF="parsnp_filtered.vcf"                 # quality-filtered VCF fed into snakemake
SAHL_TABLE="data/reference/Sahl_et_al_Table3.xlsx" # Sahl et al. SNP reference table
CLADES_TSV="config/clades.tsv"                     # clade-defining SNP file (may be updated)
SMK_FILE="snakefile_publication.smk"                # snakemake workflow
SNAKEMAKE_CORES=4                                   # CPU cores for snakemake
PARSNP_THREADS=8                                    # CPU threads for parsnp
CONDA_ENV="anthracis-pipeline"                      # conda environment name

# ---------------------------------------------------------------------------
# HELPER
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# STEP 0 - sanity checks
# ---------------------------------------------------------------------------
log "Checking required files and directories..."

[[ -d "$GENOME_DIR" ]]  || die "genome directory not found: $GENOME_DIR"
[[ -f "$REFERENCE" ]]   || die "reference genome not found: $REFERENCE"
[[ -f "$SMK_FILE" ]]    || die "snakemake file not found: $SMK_FILE"
[[ -f "$SAHL_TABLE" ]]  || die "Sahl reference table not found: $SAHL_TABLE"
[[ -f "$CLADES_TSV" ]]  || die "clades file not found: $CLADES_TSV"

FASTA_COUNT=$(find "$GENOME_DIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) | wc -l)
(( FASTA_COUNT > 0 )) || die "no genome files (.fna/.fasta/.fa) found in $GENOME_DIR"
log "Found $FASTA_COUNT genome file(s) in $GENOME_DIR."

# ---------------------------------------------------------------------------
# STEP 1 - activate the conda environment
# ---------------------------------------------------------------------------
log "Activating conda environment: $CONDA_ENV"

CONDA_BASE=$(conda info --base 2>/dev/null) \
    || die "conda not found. Please install conda first (see README.md)."

export PATH="${CONDA_BASE}/envs/${CONDA_ENV}/bin:${PATH}"

for tool in parsnp python snakemake augur; do
    command -v "$tool" &>/dev/null \
        || die "'$tool' not found. Run: bash install.sh"
done

log "Environment active — all tools found."

# ---------------------------------------------------------------------------
# STEP 2 - run parsnp on all genomes in data/genomes/
# ---------------------------------------------------------------------------
log "Running parsnp on $GENOME_DIR ..."

rm -rf "$PARSNP_OUTDIR"
rm -f "$FILTERED_VCF"

parsnp \
    -p "$PARSNP_THREADS" -v --vcf -C 1000 -e -u \
    -r "$REFERENCE" \
    -d "$GENOME_DIR" \
    -o "$PARSNP_OUTDIR" \
    -c \
    --vcf \
    --no-partition

NEW_PARSNP_VCF="${PARSNP_OUTDIR}/parsnp.vcf"
[[ -f "$NEW_PARSNP_VCF" ]] || die "parsnp did not produce $NEW_PARSNP_VCF"
log "parsnp finished. Output: $NEW_PARSNP_VCF"

# ---------------------------------------------------------------------------
# STEP 3 - filter parsnp VCF by quality (remove N / LCB variants)
# ---------------------------------------------------------------------------
log "Filtering VCF by quality..."

python3 scripts/filter_vcf_by_quality.py \
    --input  "$NEW_PARSNP_VCF" \
    --output "$FILTERED_VCF" \
    --report vcf_filter_report.txt

[[ -f "$FILTERED_VCF" ]] || die "filter_vcf_by_quality.py did not produce $FILTERED_VCF"
log "Filtered VCF written to: $FILTERED_VCF"

# ---------------------------------------------------------------------------
# STEP 4 - validate and update clade-defining SNPs against filtered VCF
# ---------------------------------------------------------------------------
log "Checking clade-defining SNPs against filtered VCF..."

set +e
python3 scripts/check_and_update_clades.py \
    --vcf    "$FILTERED_VCF" \
    --clades "$CLADES_TSV" \
    --sahl   "$SAHL_TABLE"
CLADE_EXIT=$?
set -e

if [[ $CLADE_EXIT -eq 0 ]]; then
    log "Clade SNP check passed — all defining SNPs present in VCF."
elif [[ $CLADE_EXIT -eq 2 ]]; then
    log "WARNING: some clades could not be resolved (no Sahl alternative found in VCF)."
    log "         Check config/clades_update_log.tsv for details."
    log "         The pipeline will continue. Unresolved clades will be absent from the tree."
else
    die "check_and_update_clades.py failed with unexpected exit code $CLADE_EXIT"
fi

# ---------------------------------------------------------------------------
# STEP 5 - create symlink where snakemake expects it
#           The .smk reads from 01_parsnp_out/parsnp.vcf
#           Symlink is removed at the end, leaving the legacy file untouched
# ---------------------------------------------------------------------------
log "Linking filtered VCF for snakemake..."

cp 01_parsnp_out/parsnp.vcf 01_parsnp_out/parsnp_legacy_backup.vcf
ln -sf "$(realpath "$FILTERED_VCF")" 01_parsnp_out/parsnp.vcf

# ---------------------------------------------------------------------------
# STEP 6 - run snakemake
# ---------------------------------------------------------------------------
log "Running snakemake pipeline..."

snakemake \
    --snakefile "$SMK_FILE" \
    --cores "$SNAKEMAKE_CORES" \
    --rerun-incomplete

# ---------------------------------------------------------------------------
# CLEANUP - restore legacy VCF from backup
# ---------------------------------------------------------------------------
mv 01_parsnp_out/parsnp_legacy_backup.vcf 01_parsnp_out/parsnp.vcf

log "========================================================"
log "Pipeline complete!"
log "Final output: 08_augur_out_export2auspice.json"
log "Load it in Auspice: https://auspice.us"
log "========================================================"