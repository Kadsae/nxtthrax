# NXTTHRAX - Nextstrain Phylogenetic Pipeline for *Bacillus anthracis*

Automated pipeline for phylogeographic analysis of *B. anthracis* whole-genome sequencing data. It runs **parsnp** for variant calling, merges results with a legacy VCF via **bcftools**, optionally downloads genomes from NCBI, performs quality assessment with **CheckM** and **BUSCO**, and feeds the result into an **augur / snakemake** workflow to produce an interactive **Auspice** visualisation.

---

## Features

**Automated quality control** - Assess genome completeness with CheckM and BUSCO  
**Variant calling** - Parsnp whole-genome alignment against reference (*B. anthracis* Ames-Ancestor)  
**Phylogenetic inference** - Nextstrain/Augur pipeline for tree construction and refinement  
**Clade assignment** - Automatic SNP-based clade designation  
**Interactive visualisation** - Auspice phylogeographic maps and trees  
**Flexible workflows** - Quick analysis (test set + new genomes) or full analysis (all genomes)  
**NCBI integration** - Optional genome download from NCBI  

---

## Requirements

### System
- **OS:** Linux (Ubuntu 20.04+), macOS, or Windows (WSL2)
- **RAM:** 8 GB minimum, 16 GB recommended (more for large genome sets >100 genomes)
- **Disk:** 10 GB free space minimum
- **CPU:** Multi-core recommended (4+ cores; 8+ cores for faster analysis)
- **Package manager:** [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)

### Software (installed automatically via `install.sh`)
- Python 3.11, parsnp 2.1.5, bcftools 1.23, snakemake, augur
- CheckM: Genome completeness and contamination assessment
- BUSCO: Benchmarking single-copy orthologs for completeness
- Python packages: pandas, openpyxl, vcfpy, requests, biopython
- R + ggplot2, dplyr, tidyr (for coverage plots, optional)
- Optional: NCBI datasets CLI (for genome download)

---

## Installation

### 1. Install Miniconda (if not already installed)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 2. Clone the repository

```bash
git clone https://github.com/maximilianfmayerhoferrochel/nxtthrax.git
cd nxtthrax
```

### 3. Run the installation script

```bash
bash install.sh
```

This creates a conda environment called `anthracis-pipeline` with all required tools, including CheckM and BUSCO.

**Installation time:** ~5-15 minutes (depends on internet speed and CPU)

---

## Repository Structure

```
nxtthrax/
├── install.sh                        ← Run once to set up the environment
├── run_full_analysis.sh              ← Full workflow: all genomes → tree (--with-qc runs CheckM + BUSCO)
├── run_quick_analysis.sh             ← Quick workflow: test + new genomes → tree (--with-qc runs CheckM + BUSCO)
├── environment.yml                   ← Conda environment definition (includes CheckM, BUSCO, dependencies)
├── snakefile_publication.smk         ← Snakemake/augur workflow (do not move)
├── parsnp_header.txt                 ← VCF header template (do not move)
├── exclude_samples.txt               ← Legacy sample IDs to exclude
│
├── data/
│   ├── genomes/                      ← PUT YOUR GENOME FILES HERE (.fna)
│   └── reference/
│       ├── Ames-Ancestor.fna         ← Reference genome (do not delete)
│       ├── Ames-Ancestor.gbk
│       ├── Ames-Ancestor.gff3
│       └── Sahl_et_al_Table3.xlsx    ← Clade SNP reference table (do not delete)
│
├── config/
│   ├── metadata.tsv                  ← EDIT: strain metadata, dates, locations
│   ├── lat_longs.tsv                 ← EDIT: geographic coordinates
│   ├── colors.tsv                    ← Node colour scheme for Auspice
│   └── clades.tsv                    ← Clade definitions (auto-updated by pipeline)
│
├── scripts/
│   ├── run_checkm_qc.py              ← CheckM quality assessment (do not move)
│   ├── run_busco_qc.py               ← BUSCO quality assessment (do not move)
│   ├── check_and_update_clades.py    ← Clade SNP validation (do not move)
│   ├── filter_vcf_by_quality.py      ← Quality filter for VCF
│   ├── filter_unique_snps.py         ← Removes SNPs unique to user genomes
│   ├── download_genomes.py           ← Download genomes from NCBI
│   ├── coverage_plot_ggplot.R        ← Coverage plot (R)
│   └── MST_script.r                  ← Minimum spanning tree (R)
│
├── 00_anthracis_check/               ← Local BLAST check for B. anthracis identity
│   ├── anthracis_check.sh
│   ├── dhp61/                        ← BLAST DB: chromosomal marker DHP61
│   └── pl3/                          ← BLAST DB: chromosomal marker PL3
│
└── 01_parsnp_out/
    └── parsnp.vcf                    ← Legacy VCF
```

> **Test genomes:** `data/genomes/` contains a small set of example `.fna` files with matching entries in `config/metadata.tsv`. Running `bash run_quick_analysis.sh` on these is the fastest way to verify your installation end-to-end.

---

## Quick Test

To verify your installation works end-to-end, run the quick analysis on the included test genomes:

```bash
bash run_quick_analysis.sh
```

**Expected runtime:** 10-20 minutes (without QC). With `--with-qc`: 25-50 minutes. The output `08_augur_out_export2auspice.json` can be viewed at https://auspice.us.

---

## Usage

### Workflow A: Quick Analysis (Test + New Genome → Tree)

Use this to quickly see where a new genome lands in the existing phylogeny. Parsnp runs only on the small test genome set plus your new genome(s), so it is much faster and less memory-intensive than the full analysis.

#### Basic usage (without quality control):

```bash
# 1. Place your new genome file in data/genomes/
cp your_genome.fna data/genomes/

# 2. Add a row for the new sample to config/metadata.tsv
#    Required columns: strain, date (YYYY-MM-DD), country

# 3. Run
bash run_quick_analysis.sh
```

#### With quality control (Recommended for publication):

```bash
bash run_quick_analysis.sh --with-qc
```

This runs CheckM and BUSCO quality assessments first, updates your metadata, and reports any genomes that fail quality thresholds.

**Steps performed:**

| Step | Tool | What happens |
|------|------|--------------|
| 0 (Optional) | CheckM + BUSCO | Assesses genome completeness/contamination (CheckM) and conserved orthologs (BUSCO) |
| 1 | parsnp | Aligns test genomes + your new genome(s) against reference → `01_parsnp_out_new/parsnp.vcf` |
| 2 | bcftools merge | Merges new VCF with legacy VCF, excluding samples in `exclude_samples.txt` |
| 3 | filter_vcf_by_quality.py | Removes variants with FILTER starting with `N` (N, N:LCB, N:ALN, ...) |
| 4 | filter_unique_snps.py | Removes SNPs present only in the user's new genome(s) |
| 5 | check_and_update_clades.py | Validates clade-defining SNPs; auto-updates `config/clades.tsv` if needed |
| 6 | snakemake / augur | Full phylogenetic workflow: tree → refine → ancestral → traits → clades → export |

**Total runtime:**
- Without QC: ~10-20 minutes
- With CheckM only: ~15-30 minutes
- With BUSCO only: ~20-50 minutes
- With both CheckM + BUSCO: ~30-80 minutes
(Times depend on genome count, size, and CPU cores. Use `--checkm-threads` and `--busco-threads` to parallelize.)

---

### Workflow B: Full Analysis (All Genomes → Tree)

Use this when you want to rebuild the complete phylogeny with all your genomes in `data/genomes/`. Parsnp aligns everything from scratch, which requires more time and RAM.

#### Basic usage:

```bash
bash run_full_analysis.sh
```

#### With quality control (Recommended for publication):

```bash
bash run_full_analysis.sh --with-qc
```

#### Advanced options:

```bash
# CheckM only (skip BUSCO)
bash run_full_analysis.sh --with-qc --skip-busco

# BUSCO only (skip CheckM)
bash run_full_analysis.sh --with-qc --skip-checkm

# Strict quality control with both tools
bash run_full_analysis.sh --with-qc \
    --checkm-min-completeness 95 \
    --checkm-max-contamination 5 \
    --busco-min-completeness 95

# Use specific BUSCO lineage (for non-B.anthracis genomes)
bash run_full_analysis.sh --with-qc \
    --busco-lineage bacteria_odb10

# Faster QC with more cores
bash run_full_analysis.sh --with-qc \
    --checkm-threads 16 \
    --busco-threads 16
```

**Total runtime:**
- Without QC: ~20-60 minutes (depends on genome count)
- With CheckM only: ~25-90 minutes
- With BUSCO only: ~40-120 minutes
- With both CheckM + BUSCO: ~60-180 minutes

**Command line options:**

```
--with-qc                           Enable CheckM + BUSCO quality control (default: off)
--skip-checkm                       Skip CheckM (use with --with-qc to run BUSCO only)
--skip-busco                        Skip BUSCO (use with --with-qc to run CheckM only)
--checkm-min-completeness PCT       Minimum genome completeness (default: 90%)
--checkm-max-contamination PCT      Maximum contamination (default: 10%)
--checkm-threads N                  CPU threads for CheckM (default: 4)
--busco-lineage LINEAGE             BUSCO lineage dataset (default: bacillales_odb10)
--busco-min-completeness PCT        Minimum BUSCO completeness (default: 90%)
--busco-threads N                   CPU threads for BUSCO (default: 4)
--help                              Show usage information
```

---

## Quality Control Assessment

### Tools

Your pipeline includes two complementary quality control tools:

1. **CheckM** - Assesses genome completeness and contamination by identifying marker genes
2. **BUSCO** - Assesses genome completeness by searching for conserved single-copy orthologs

### What is CheckM?

CheckM estimates genome quality by assessing **completeness** (percentage of expected genes present) and **contamination** (percentage of unexpected extra genes). It's essential for publication-quality phylogenomic work.

### What is BUSCO?

BUSCO (Benchmarking Universal Single-Copy Orthologs) evaluates genome completeness by searching for lineage-specific conserved single-copy genes. Unlike CheckM (which uses marker genes), BUSCO looks for complete gene families expected to be present in a genome.

For *B. anthracis*, the default lineage is `bacillales_odb10` (specific to Bacillales order), which is faster and more specific than the broader `bacteria_odb10` dataset.

### Run Both Tools (Recommended for Publication)

```bash
bash run_full_analysis.sh --with-qc
```

This runs **both** CheckM and BUSCO, providing complementary quality metrics.

### Run Individual Tools

```bash
# CheckM only
bash run_full_analysis.sh --with-qc --skip-busco

# BUSCO only
bash run_full_analysis.sh --with-qc --skip-checkm
```

### BUSCO Lineage Options

For *B. anthracis* (default and recommended):
```bash
--busco-lineage bacillales_odb10
```

For other gram-positive bacteria:
```bash
--busco-lineage bacteria_odb10
```

For automatic detection (slower, no download needed):
```bash
--busco-lineage auto-lineage-prok
```

### Quick Quality Assessment

Run CheckM independently:

```bash
python scripts/run_checkm_qc.py \
    --genome-dir data/genomes \
    --metadata config/metadata.tsv \
    --update-metadata
```

Run BUSCO independently:

```bash
python scripts/run_busco_qc.py \
    --genome-dir data/genomes \
    --lineage bacillales_odb10 \
    --metadata config/metadata.tsv \
    --update-metadata
```

**Output files (CheckM):**
- `checkm_quality_report.tsv` - Summary table of all genomes
- `genomes_failed_qc.tsv` - List of genomes failing quality thresholds (if any)
- `checkm_results/` - Full CheckM output directory

**Output files (BUSCO):**
- `busco_quality_report.tsv` - Summary table of all genomes
- `genomes_failed_busco.tsv` - List of genomes failing quality thresholds (if any)
- `busco_results/` - Full BUSCO output directory per genome

### Recommended Quality Thresholds

For peer-reviewed publication, use these thresholds:

| Tool | Threshold | Completeness | Contamination | Use Case |
|------|-----------|--------------|---------------|----------|
| CheckM | Strict | ≥95% | ≤5% | High-quality isolates only |
| CheckM | **Standard** | **≥90%** | **≤10%** | **Recommended for publication** |
| CheckM | Permissive | ≥80% | ≤15% | Include borderline genomes |
| BUSCO | Strict | ≥95% | N/A | High-quality isolates only |
| BUSCO | **Standard** | **≥90%** | **N/A** | **Recommended for publication** |
| BUSCO | Permissive | ≥80% | N/A | Include borderline genomes |

**Recommended:** Use both CheckM (90%/≤10%) and BUSCO (≥90%) thresholds for publication-quality work.

### Interpreting Results

From `checkm_quality_report.tsv`:

```
genome            completeness  contamination
sample_1.fna      99.5          0.3          ✓ Excellent
sample_2.fna      92.1          4.2          ✓ Good
sample_3.fna      78.2          12.1         ✗ Failed (below 90% completeness)
```

From `busco_quality_report.tsv`:

```
genome            completeness  complete_single_copy  fragmented  missing
sample_1.fna      98.7          320                   1           2       ✓ Excellent
sample_2.fna      91.3          296                   3           8       ✓ Good
sample_3.fna      78.5          254                   12          47      ✗ Failed (below 90% completeness)
```

---

### Optional: Download Genomes from NCBI

If your `metadata.tsv` contains NCBI accession numbers (GCF_/GCA_ format):

```bash
# Requires: conda install -c conda-forge ncbi-datasets-cli
python scripts/download_genomes.py \
    --metadata config/metadata.tsv \
    --outdir data/genomes \
    --retry 3
```

---

### B. anthracis Identity Check (Optional)

Before adding a new genome, confirm it is *B. anthracis* by running a local BLAST against two chromosomal markers (DHP61 and PL3):

```bash
cd 00_anthracis_check
# Edit anthracis_check.sh to point to your genome directory, then:
bash anthracis_check.sh
```

Results are written to `dhp61_results.txt` and `pl3_results.txt`.  
See [NCBI BLAST setup instructions](https://www.ncbi.nlm.nih.gov/books/NBK279690/) if you need to configure a local BLAST database.

---

### Coverage & Minimum Spanning Tree Plots (Optional)

```bash
# Coverage plot (requires R)
Rscript scripts/coverage_plot_ggplot.R

# Minimum spanning tree
Rscript scripts/MST_script.r
# Upload resulting tree file to https://achtman-lab.github.io/GrapeTree/MSTree_holder.html
```

---

## Output

### Main Output

| File | Description |
|------|-------------|
| `08_augur_out_export2auspice.json` | **Main output** - phylogenetic tree with metadata, ready for Auspice |
| `02_augur_out_tree_raw.nwk` | Raw phylogenetic tree (Newick format) |
| `03_augur_out_tree_refined.nwk` | Time-refined phylogenetic tree (Newick format) |

### Quality Control Output (if --with-qc used)

| File | Description |
|------|-------------|
| `checkm_quality_report.tsv` | CheckM summary table of completeness and contamination |
| `genomes_failed_qc.tsv` | CheckM: list of genomes failing thresholds (optional) |
| `checkm_results/` | CheckM full output and marker gene analysis |
| `busco_quality_report.tsv` | BUSCO summary table of completeness |
| `genomes_failed_busco.tsv` | BUSCO: list of genomes failing thresholds (optional) |
| `busco_results/` | BUSCO full output per genome |

### VCF and Filtering Output

| File | Description |
|------|-------------|
| `parsnp_new_merged.vcf` | Merged VCF (new + legacy) |
| `parsnp_filtered.vcf` | Quality-filtered VCF |
| `parsnp_filtered_unique.vcf` | Unique-SNP-filtered VCF (quick analysis only) |
| `vcf_filter_report.txt` | VCF filtering statistics |

### Clade and Metadata Output

| File | Description |
|------|-------------|
| `config/clades_update_log.tsv` | Log of any clade SNP changes |
| `config/metadata.tsv` | Updated with CheckM metrics (if --update-metadata used) |

### Visualisation

To view your results in **Auspice**:

1. **Online:** Go to https://auspice.us and drag-and-drop `08_augur_out_export2auspice.json`
2. **Locally:** 
```bash
nextstrain view auspice/ --port 4001
# Open http://localhost:4001 in your browser
```

---

## Configuration Files

### `config/metadata.tsv`

Tab-separated strain metadata with quality control metrics.

**Required columns:**
- `strain` - Unique genome identifier (matches filename without .fna)
- `date` - Collection date in YYYY-MM-DD format (or YYYY-XX-XX if unknown)
- `country` - Country of origin

**Optional columns:**
- `region` - Region/province/state within country
- `host` - Host species
- `isolation_source` - Source material (e.g., clinical, environmental)
- `checkm_completeness` - CheckM completeness (%) - auto-populated if --update-metadata used
- `checkm_contamination` - CheckM contamination (%) - auto-populated if --update-metadata used
- `checkm_strain_heterogeneity` - CheckM strain heterogeneity - auto-populated if --update-metadata used
- `busco_completeness` - BUSCO completeness (%) - auto-populated if --update-metadata used
- `busco_complete_single` - BUSCO single-copy complete count - auto-populated if --update-metadata used
- `busco_complete_duplicated` - BUSCO duplicated complete count - auto-populated if --update-metadata used
- `busco_fragmented` - BUSCO fragmented count - auto-populated if --update-metadata used
- `busco_missing` - BUSCO missing count - auto-populated if --update-metadata used

Example:
```tsv
strain              date            country     checkm_completeness    checkm_contamination    busco_completeness
sample_001.fna      2020-05-15      Germany     99.52                  0.29                    98.7
sample_002.fna      2019-03-22      France      92.18                  4.21                    91.3
```

### `config/lat_longs.tsv`

Geographic coordinates for map visualisation in Auspice.

### `config/colors.tsv`

Custom colour scheme for Auspice nodes (optional).

### `config/clades.tsv`

Clade definitions based on SNP patterns. Auto-updated by the pipeline using SNP references from `Sahl_et_al_Table3.xlsx`.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `conda: command not found` | Install Miniconda, then `source ~/.bashrc` |
| `parsnp: command not found` | Re-run `bash install.sh` or activate env: `conda activate anthracis-pipeline` |
| `checkm: command not found` | Run `bash install.sh` again or `conda install -c bioconda checkm-genome` |
| `vcfpy not found` | `conda activate anthracis-pipeline && pip install vcfpy` |
| `Permission denied` | `chmod +x install.sh run_full_analysis.sh run_quick_analysis.sh` |
| Parsnp is very slow | Expected - whole-genome alignment is compute-intensive. Use more threads: edit `PARSNP_THREADS` in the script or wait. |
| CheckM is very slow | Expected - marker gene identification takes time. Use more threads: `--checkm-threads 16` |
| BUSCO is very slow | Expected - ortholog searching takes time. Use `--busco-threads 16` or a faster lineage like `bacteria_odb10` |
| BUSCO not found | Run `bash install.sh` again or `conda install -c bioconda busco` |
| BUSCO failed / skipped | Install BUSCO (`conda install -c bioconda busco`) and retry with `--with-qc` |
| Some clades not resolved | Non-fatal warning. Check `config/clades_update_log.tsv`. Manually add alternative SNPs to `config/clades.tsv` if needed. |
| Out of memory | Try a smaller subset first. ~50 genomes require 8 GB RAM; 500 genomes require 16+ GB RAM. |
| Results differ between runs | Normal - tree inference uses randomisation. Trees should be topologically similar. |
| CheckM failed / skipped | Install CheckM (`conda install -c bioconda checkm-genome`) and retry with `--with-qc` |

---

## References

- **Nextstrain / Augur / Auspice:** https://nextstrain.org
- **ParSNP v2.0:** https://github.com/marbl/parsnp
- **BCFtools:** https://samtools.github.io/bcftools/
- **Snakemake:** https://snakemake.readthedocs.io/
- **CheckM:** https://github.com/ecogenomics/checkm
- **BUSCO:** https://busco.ezlab.org/
- **Sahl et al. (2016):** Clade SNP reference

---

## License

See [LICENSE](LICENSE) for details.

---

## Support & Feedback

If you use NXTTHRAX in your research, please cite us!

For issues, questions, or feature requests, please open an issue on GitHub or contact the authors.

Last updated: June 2026