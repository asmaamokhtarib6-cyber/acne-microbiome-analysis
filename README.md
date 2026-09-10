# Acne Skin Microbiome Analysis

A learning project where I worked through a full 16S rRNA amplicon
sequencing pipeline, from raw public sequencing data to statistics in R.

## About

I built this to actually understand microbiome bioinformatics, not just
read about it. It reprocesses public data from Sun et al. (2025) as
practice, so this is a learning exercise, not a research paper. No new
scientific claims here, just me working through a real pipeline
end to end.

## What I learned

- QIIME2: importing raw reads, running DADA2 to denoise them, then
  assigning taxonomy
- Pulling and handling public sequencing data from NCBI SRA using
  Entrez Direct
- R statistics: Wilcoxon tests, PERMANOVA, alpha and beta diversity,
  Benjamini-Hochberg correction, hierarchical clustering
- Working on a remote Linux server over SSH, and keeping long jobs
  alive with `screen`

## Pipeline

1. Download raw samples (FASTQ) from NCBI SRA
2. Import into QIIME2
3. Check read quality
4. Run DADA2 (denoise, build the feature table)
5. Assign taxonomy against the SILVA database
6. Statistical analysis in R

## Data source

- Raw data: NCBI SRA, BioProject [PRJNA995545](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA995545)
- Original paper: Sun et al. 2025, *Microbiology Spectrum*,
  DOI: [10.1128/spectrum.02603-24](https://doi.org/10.1128/spectrum.02603-24)

## Repository structure

- `scripts/` bash scripts for downloading from SRA and running the QIIME2 pipeline
- `analysis/` the R script for statistical analysis
- `figures/` output plots: PCoA, heatmap, pipeline diagrams

## How to run

1. `scripts/01_download_sra.sh` downloads raw FASTQ files from NCBI
2. `scripts/02_qiime2_pipeline.sh` runs QIIME2 import, DADA2, taxonomy
3. `analysis/sun2025_analysis.R` runs the statistical analysis (R or Google Colab)

## Note

This repository doesn't include data from a related project that uses
unpublished patient data. That data stays private and isn't part of
this repo.

## Author

Asmaa Mokhtar, MSc Biomedical Science, Bioinformatics track
