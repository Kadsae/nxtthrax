#!/bin/bash

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
SMK_FILE="snakefile_publication.smk"                # snakemake workflow (do not move)
SNAKEMAKE_CORES=4                                   # CPU cores for snakemake
PARSNP_THREADS=8                                    # CPU threads for parsnp
CONDA_ENV="anthracis-pipeline"                      # conda environment name

# QC / CheckM settings
CHECKM_ENABLED=false                                # default: disabled
CHECKM_MIN_COMPLETENESS=90                          # minimum genome completeness (%)
CHECKM_MAX_CONTAMINATION=10                         # maximum contamination (%)
CHECKM_THREADS=4                                    # CPU threads for CheckM

# QC / BUSCO settings
BUSCO_ENABLED=false                                 # default: disabled
BUSCO_LINEAGE="bacillales_odb10"                    # BUSCO lineage for B. anthracis
BUSCO_MIN_COMPLETENESS=90                           # minimum BUSCO completeness (%)
BUSCO_THREADS=4                                     # CPU threads for BUSCO

# QC Tool skip flags (when --with-qc is used)
SKIP_CHECKM=false                                   # skip CheckM even if --with-qc
SKIP_BUSCO=false                                    # skip BUSCO even if --with-qc

# ---------------------------------------------------------------------------
# HELPER
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# PARSE COMMAND LINE ARGUMENTS
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --with-qc)
            CHECKM_ENABLED=true
            BUSCO_ENABLED=true
            shift
            ;;
        --skip-checkm)
            SKIP_CHECKM=true
            shift
            ;;
        --skip-busco)
            SKIP_BUSCO=true
            shift
            ;;
        --checkm-min-completeness)
            CHECKM_MIN_COMPLETENESS="$2"
            shift 2
            ;;
        --checkm-max-contamination)
            CHECKM_MAX_CONTAMINATION="$2"
            shift 2
            ;;
        --checkm-threads)
            CHECKM_THREADS="$2"
            shift 2
            ;;
        --busco-lineage)
            BUSCO_LINEAGE="$2"
            shift 2
            ;;
        --busco-min-completeness)
            BUSCO_MIN_COMPLETENESS="$2"
            shift 2
            ;;
        --busco-threads)
            BUSCO_THREADS="$2"
            shift 2
            ;;
        --help)
            cat << 'HELP'
Full phylogenetic analysis workflow for NXTTHRAX

Usage:
    bash run_full_analysis.sh [OPTIONS]

Options:
    --with-qc                           Enable CheckM + BUSCO quality control (default: off)
    --skip-checkm                       Skip CheckM (use with --with-qc to run BUSCO only)
    --skip-busco                        Skip BUSCO (use with --with-qc to run CheckM only)
    --checkm-min-completeness PCT       Minimum genome completeness (default: 90)
    --checkm-max-contamination PCT      Maximum contamination (default: 10)
    --checkm-threads N                  CPU threads for CheckM (default: 4)
    --busco-lineage LINEAGE             BUSCO lineage dataset (default: bacillales_odb10)
    --busco-min-completeness PCT        Minimum BUSCO completeness (default: 90)
    --busco-threads N                   CPU threads for BUSCO (default: 4)
    --help                              Show this help message

Examples:
    # Standard run (no QC)
    bash run_full_analysis.sh

    # With both CheckM and BUSCO
    bash run_full_analysis.sh --with-qc

    # CheckM only (skip BUSCO)
    bash run_full_analysis.sh --with-qc --skip-busco

    # BUSCO only (skip CheckM)
    bash run_full_analysis.sh --with-qc --skip-checkm

    # Strict quality control for both tools
    bash run_full_analysis.sh --with-qc \
        --checkm-min-completeness 95 \
        --checkm-max-contamination 5 \
        --busco-min-completeness 95

    # Use more cores for faster analysis
    bash run_full_analysis.sh --with-qc --checkm-threads 16 --busco-threads 16
HELP
            exit 0
            ;;
        *)
            die "Unknown option: $1. Use --help for usage information."
            ;;
    esac
done

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
# STEP 1.5 - Optional: Run quality assessment (CheckM and/or BUSCO)
# ---------------------------------------------------------------------------
if [[ "$CHECKM_ENABLED" == "true" ]] || [[ "$BUSCO_ENABLED" == "true" ]]; then
    log "=========================================================="
    log "STEP 1.5: Quality Control Assessment (Enabled)"
    log "=========================================================="
    
    # Adjust tool flags based on skip settings
    if [[ "$SKIP_CHECKM" == "true" ]]; then
        CHECKM_ENABLED=false
    fi
    if [[ "$SKIP_BUSCO" == "true" ]]; then
        BUSCO_ENABLED=false
    fi
    
    # Run CheckM if enabled
    if [[ "$CHECKM_ENABLED" == "true" ]]; then
        log "Running CheckM on genomes (this may take 5-30 minutes)..."
        
        if command -v checkm &>/dev/null; then
            python3 scripts/run_checkm_qc.py \
                --genome-dir "$GENOME_DIR" \
                --output-dir checkm_results \
                --metadata config/metadata.tsv \
                --min-completeness "$CHECKM_MIN_COMPLETENESS" \
                --max-contamination "$CHECKM_MAX_CONTAMINATION" \
                --update-metadata \
                --threads "$CHECKM_THREADS" \
                || log "WARNING: CheckM failed. Proceeding with analysis."
            
            # Check if any genomes failed QC
            if [[ -f genomes_failed_qc.tsv ]]; then
                failed_count=$(tail -n +2 genomes_failed_qc.tsv 2>/dev/null | wc -l)
                if (( failed_count > 0 )); then
                    log "WARNING: $failed_count genome(s) failed CheckM thresholds"
                    log "See genomes_failed_qc.tsv for details"
                    log "Consider removing these genomes from $GENOME_DIR before proceeding"
                fi
            fi
            
            log "CheckM quality control complete ✓"
        else
            log "WARNING: CheckM not found. Install with: conda install -c bioconda checkm-genome"
            log "Skipping CheckM and proceeding with analysis."
        fi
    fi
    
    # Run BUSCO if enabled
    if [[ "$BUSCO_ENABLED" == "true" ]]; then
        log "Running BUSCO on genomes (this may take 10-60 minutes)..."
        
        if command -v busco &>/dev/null; then
            python3 scripts/run_busco_qc.py \
                --genome-dir "$GENOME_DIR" \
                --output-dir busco_results \
                --metadata config/metadata.tsv \
                --lineage "$BUSCO_LINEAGE" \
                --min-completeness "$BUSCO_MIN_COMPLETENESS" \
                --update-metadata \
                --threads "$BUSCO_THREADS" \
                || log "WARNING: BUSCO failed. Proceeding with analysis."
            
            # Check if any genomes failed QC
            if [[ -f genomes_failed_busco.tsv ]]; then
                failed_count=$(tail -n +2 genomes_failed_busco.tsv 2>/dev/null | wc -l)
                if (( failed_count > 0 )); then
                    log "WARNING: $failed_count genome(s) failed BUSCO thresholds"
                    log "See genomes_failed_busco.tsv for details"
                    log "Consider removing these genomes from $GENOME_DIR before proceeding"
                fi
            fi
            
            log "BUSCO quality control complete ✓"
        else
            log "WARNING: BUSCO not found. Install with: conda install -c bioconda busco"
            log "Skipping BUSCO and proceeding with analysis."
        fi
    fi
fi

# ---------------------------------------------------------------------------
# STEP 2 - run parsnp on all genomes in data/genomes/
# ---------------------------------------------------------------------------
log "=========================================================="
log "STEP 2: Running ParSNP"
log "=========================================================="
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
log "parsnp finished. Output: $NEW_PARSNP_VCF ✓"

# ---------------------------------------------------------------------------
# STEP 3 - filter parsnp VCF by quality (remove N / LCB variants)
# ---------------------------------------------------------------------------
log "=========================================================="
log "STEP 3: VCF Quality Filtering"
log "=========================================================="
log "Filtering VCF by quality..."

python3 scripts/filter_vcf_by_quality.py \
    --input  "$NEW_PARSNP_VCF" \
    --output "$FILTERED_VCF" \
    --report vcf_filter_report.txt

[[ -f "$FILTERED_VCF" ]] || die "filter_vcf_by_quality.py did not produce $FILTERED_VCF"
log "Filtered VCF written to: $FILTERED_VCF ✓"

# ---------------------------------------------------------------------------
# STEP 4 - validate and update clade-defining SNPs against filtered VCF
# ---------------------------------------------------------------------------
log "=========================================================="
log "STEP 4: Clade SNP Validation"
log "=========================================================="
log "Checking clade-defining SNPs against filtered VCF..."

set +e
python3 scripts/check_and_update_clades.py \
    --vcf    "$FILTERED_VCF" \
    --clades "$CLADES_TSV" \
    --sahl   "$SAHL_TABLE"
CLADE_EXIT=$?
set -e

if [[ $CLADE_EXIT -eq 0 ]]; then
    log "Clade SNP check passed — all defining SNPs present in VCF. ✓"
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
log "=========================================================="
log "STEP 5: Phylogenetic Analysis (Snakemake/Augur)"
log "=========================================================="
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

log "=========================================================="
log "PIPELINE COMPLETE!"
log "=========================================================="
log "Final output: 08_augur_out_export2auspice.json"
log "Load it in Auspice: https://auspice.us"

if [[ "$CHECKM_ENABLED" == "true" ]] || [[ "$BUSCO_ENABLED" == "true" ]]; then
    log ""
    log "Quality Control Summary:"
    if [[ "$CHECKM_ENABLED" == "true" ]]; then
        log "  - CheckM results: checkm_results/"
        log "  - CheckM report: checkm_quality_report.tsv"
    fi
    if [[ "$BUSCO_ENABLED" == "true" ]]; then
        log "  - BUSCO results: busco_results/"
        log "  - BUSCO report: busco_quality_report.tsv"
    fi
    log "  - Metadata updated: config/metadata.tsv"
fi

log "=========================================================="
