#!/usr/bin/env bash

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
ENV_NAME="anthracis-pipeline"
ENV_FILE="environment.yml"

# ---------------------------------------------------------------------------
# HELPER
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# STEP 0 - find conda or mamba
# ---------------------------------------------------------------------------
log "Looking for conda or mamba..."

if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    die "conda or mamba not found. Install Miniconda first:
  https://docs.conda.io/en/latest/miniconda.html
  Then restart your terminal and run this script again."
fi

log "Using: $CONDA_CMD ($($CONDA_CMD --version))"

# ---------------------------------------------------------------------------
# STEP 1 - create or update the conda environment
# ---------------------------------------------------------------------------
log "Setting up environment '$ENV_NAME' from $ENV_FILE..."

[[ -f "$ENV_FILE" ]] || die "$ENV_FILE not found. Make sure you are in the project directory."

if $CONDA_CMD env list | grep -q "^${ENV_NAME} "; then
    log "Environment '$ENV_NAME' already exists."
    read -p "Update it? (y/n) " -n 1 -r; echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        $CONDA_CMD env update -n "$ENV_NAME" --file "$ENV_FILE" --prune
        log "Environment updated."
    else
        log "Skipping environment update."
    fi
else
    log "Creating environment (this may take 5-15 minutes)..."
    $CONDA_CMD env create -n "$ENV_NAME" --file "$ENV_FILE" \
        || die "Failed to create environment."
    log "Environment created."
fi

# ---------------------------------------------------------------------------
# STEP 2 - verify tools inside the environment
# ---------------------------------------------------------------------------
log "Verifying installed tools..."

CONDA_BASE=$($CONDA_CMD info --base)
ACTIVATE_SCRIPT=$(mktemp)

cat > "$ACTIVATE_SCRIPT" << 'INNER'
#!/usr/bin/env bash
set -e
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate anthracis-pipeline

failed=0
for tool in python parsnp bcftools snakemake augur; do
    if command -v "$tool" &> /dev/null; then
        echo "  [ok] $tool"
    else
        echo "  [!!] $tool — NOT FOUND"
        failed=$((failed + 1))
    fi
done

exit $failed
INNER

chmod +x "$ACTIVATE_SCRIPT"

if bash "$ACTIVATE_SCRIPT"; then
    log "All tools verified."
else
    rm "$ACTIVATE_SCRIPT"
    die "Some tools are missing. Try re-running this script or check environment.yml."
fi

rm "$ACTIVATE_SCRIPT"

# ---------------------------------------------------------------------------
# STEP 3 - check project structure
# ---------------------------------------------------------------------------
log "Checking project structure..."

missing=0

for f in snakefile_publication.smk parsnp_header.txt exclude_samples.txt environment.yml; do
    [[ -f "$f" ]] || { log "  WARNING: missing file: $f"; missing=$((missing + 1)); }
done

for d in data/genomes data/reference config scripts; do
    [[ -d "$d" ]] || { log "  WARNING: missing directory: $d"; missing=$((missing + 1)); }
done

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------
log "========================================================"
log "Installation complete!"
if (( missing > 0 )); then
    log "WARNING: $missing file(s)/directory(ies) missing — check structure above."
fi
log "Run the pipeline:"
log "  bash run_quick_analysis.sh   (quick placement with test genomes)"
log "  bash run_full_analysis.sh    (full rebuild from all genomes)"
log "========================================================"