# NXTTHRAX — Nextstrain Phylogenetic Pipeline for *Bacillus anthracis*

Automated pipeline for phylogeographic analysis of *B. anthracis* whole-genome sequencing data. It runs **parsnp** for variant calling, merges results with a legacy VCF via **bcftools**, optionally downloads genomes from NCBI, and feeds the result into an **augur / snakemake** workflow to produce an interactive **Auspice** visualisation.

---

## Requirements

### System
- **OS:** Linux (Ubuntu 20.04+), macOS, or Windows (WSL2)
- **RAM:** 8 GB minimum, 16 GB recommended
- **Disk:** 10 GB free space
- **CPU:** Multi-core recommended (4+ cores)
- **Package manager:** [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)

### Software (installed automatically via `install.sh`)
- Python 3.11, parsnp 2.1.5, bcftools 1.23, snakemake, augur
- Python packages: pandas, openpyxl, vcfpy, requests, biopython
- R + ggplot2, dplyr, tidyr (for coverage plots)
- Optional: NCBI datasets CLI (for genome download)

---

## Installation

```bash
# 1. Install Miniconda if not already present
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# 2. Clone the repository
git clone https://github.com/maximilianfmayerhoferrochel/nxtthrax.git
cd nxtthrax

# 3. Run the installation script (creates the 'anthracis-pipeline' conda env)
bash install.sh
```

Installation takes roughly 5–15 minutes depending on internet speed.

---

## Repository Structure

```
nxtthrax/
├── install.sh                        ← Run once to set up the environment
├── run_full_analysis.sh              ← Full workflow: new genomes → tree
├── run_quick_analysis.sh             ← Quick workflow: existing VCF → tree
├── environment.yml                   ← Conda environment definition
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
    └── parsnp.vcf                    ← Legacy VCF (do not delete)
```

> **Test genomes:** `data/genomes/` contains a small set of example `.fna` files with matching entries in `config/metadata.tsv`. Running `bash run_quick_analysis.sh` on these is the fastest way to verify your installation end-to-end.

---

## Quick Test

To verify your installation works end-to-end, run the quick analysis on the included test genomes:

```bash
bash run_quick_analysis.sh
```

The test genomes in `data/genomes/` already have matching rows in `config/metadata.tsv`, so no editing is needed. Expected runtime: 10–20 minutes. The output `08_augur_out_export2auspice.json` can be viewed at https://auspice.us.

---

## Usage

### Workflow A — Quick Analysis (test genomes + new genome → tree)

Use this to quickly see where a new genome lands in the phylogeny. Parsnp runs only on the small test genome set plus your new genome, so it is much faster and less memory-intensive than the full analysis.

```bash
# 1. Place your new genome file in data/genomes/
cp your_genome.fna data/genomes/

# 2. Add a row for the new sample to config/metadata.tsv
#    Required columns: strain, date (YYYY-MM-DD), country

# 3. Run
bash run_quick_analysis.sh
```

**Steps performed:**

| Step | Tool | What happens |
|------|------|--------------|
| 1 | parsnp | Aligns test genomes + your new genome against the reference → `01_parsnp_out_new/parsnp.vcf` |
| 2 | bcftools merge | Merges new VCF with legacy VCF, excluding samples in `exclude_samples.txt` |
| 3 | filter_vcf_by_quality.py | Removes variants with FILTER starting with `N` (N, N:LCB, N:ALN, ...) |
| 4 | filter_unique_snps.py | Removes SNPs present only in the user's new genome(s) |
| 5 | check_and_update_clades.py | Validates clade-defining SNPs; auto-updates `config/clades.tsv` if needed |
| 6 | snakemake / augur | Full phylogenetic workflow: tree → refine → ancestral → translate → traits → clades → export |

**Total runtime:** ~10–20 minutes

---

### Workflow B — Full Analysis (all genomes → tree)

Use this when you want to rebuild the complete phylogeny with all your genomes in `data/genomes/`. Parsnp aligns everything from scratch, which requires more time and RAM.

```bash
bash run_full_analysis.sh
```

**Total runtime:** ~20–60 minutes (scales with genome count and available hardware)

---

### Download Genomes from NCBI (Optional)

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

| File | Description |
|------|-------------|
| `08_augur_out_export2auspice.json` | **Main output** — upload to https://auspice.us |
| `parsnp_new_merged.vcf` | Merged VCF (new + legacy) |
| `parsnp_filtered.vcf` | Quality-filtered VCF |
| `parsnp_filtered_unique.vcf` | Unique-SNP-filtered VCF (input to snakemake, quick analysis only) |
| `vcf_filter_report.txt` | Filtering statistics |
| `config/clades_update_log.tsv` | Log of any clade SNP swaps |
| `02_augur_out_tree_raw.nwk` | Raw phylogenetic tree |
| `03_augur_out_tree_refined.nwk` | Time-refined tree |

To visualise the result, go to **https://auspice.us** and drag-and-drop `08_augur_out_export2auspice.json` onto the page. No account or installation needed.

To view locally:
```bash
nextstrain view auspice/ --port 4001
# Open http://localhost:4001
```

---

## Configuration Files

### `config/metadata.tsv`
Tab-separated strain metadata. Required columns: `strain`, `date` (YYYY-MM-DD), `country`. Optional: `region`, `accession`.

### `config/lat_longs.tsv`
Geographic coordinates used for map visualisation in Auspice.

### `config/colors.tsv`
Custom colour scheme for Auspice nodes. Optional.

### `config/clades.tsv`
Clade definitions based on SNP patterns. Auto-updated by `check_and_update_clades.py` using SNP references from `Sahl_et_al_Table3.xlsx`.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `conda: command not found` | Install Miniconda, then `source ~/.bashrc` |
| `parsnp: command not found` | Re-run `bash install.sh` or activate env: `conda activate anthracis-pipeline` |
| `vcfpy not found` | `conda activate anthracis-pipeline && pip install vcfpy` |
| `Permission denied` | `chmod +x install.sh run_full_analysis.sh run_quick_analysis.sh` |
| Parsnp is slow | Expected — whole-genome alignment is compute-intensive. Use more threads: `PARSNP_THREADS=16` in `run_full_analysis.sh` |
| Some clades not resolved | Non-fatal warning. Check `config/clades_update_log.tsv`. Manually add alternative SNPs to `config/clades.tsv` if needed |
| Out of memory | Try a smaller subset first. 50 genomes ≈ 8 GB RAM; 500 genomes ≈ 16+ GB RAM |
| Results differ between runs | Normal — tree inference uses randomisation. Trees should be topologically similar |

---

## References

This workflow follows the [Nextstrain phylogenetic workflow tutorials](https://docs.nextstrain.org/en/latest/tutorials/).

- **Nextstrain / Augur / Auspice:** https://nextstrain.org
- **parsnp 2.0:** https://github.com/marbl/parsnp
- **BCFtools:** https://samtools.github.io/bcftools/
- **Snakemake:** https://snakemake.readthedocs.io/
- **Clade SNP reference:** Sahl et al. (2016)

---

## License

See [LICENSE](LICENSE) for details.