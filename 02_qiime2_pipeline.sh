#!/bin/bash
# =============================================================================
# 02_qiime2_pipeline.sh
# QIIME2 processing pipeline: import -> quality check -> DADA2 denoising ->
# taxonomic classification, for the Sun et al. 2025 dataset.
#
# Requirements: QIIME2 (tested with qiime2-amplicon-2024.10), a sample
# metadata file with columns matching your study design, and the raw_fastq/
# directory produced by 01_download_sra.sh.
#
# Usage (run from inside the QIIME2 conda environment):
#   conda activate qiime2-amplicon-2024.10
#   chmod +x 02_qiime2_pipeline.sh
#   ./02_qiime2_pipeline.sh
# =============================================================================

set -e

cd ~/sun2025_analysis

# -----------------------------------------------------------------------------
# STEP 1: Build the QIIME2 manifest file (sample-id -> FASTQ file paths).
# Expects a metadata file named sun_metadata_final.csv with at least a
# sample-id column matching the SRA accessions downloaded in step 01.
# -----------------------------------------------------------------------------
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
for f in raw_fastq/*_1.fastq.gz; do
  acc=$(basename "$f" _1.fastq.gz)
  echo -e "${acc}\t$(pwd)/raw_fastq/${acc}_1.fastq.gz\t$(pwd)/raw_fastq/${acc}_2.fastq.gz" >> manifest.tsv
done
echo "Manifest rows (expect samples + 1 header):"
wc -l manifest.tsv

# Convert metadata CSV to QIIME2's required TSV format (first column renamed
# to "sample-id"). Adjust column numbers/names to match your own metadata.
awk -F',' 'BEGIN{OFS="\t"} NR==1{$1="sample-id"} {print $1,$2,$3,$4,$5,$6,$7}' \
  sun_metadata_final.csv > sample-metadata.tsv

# -----------------------------------------------------------------------------
# STEP 2: Import into QIIME2
# -----------------------------------------------------------------------------
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# -----------------------------------------------------------------------------
# STEP 3: Summarize read quality (inspect demux.qzv at view.qiime2.org to
# choose truncation lengths for DADA2 below)
# -----------------------------------------------------------------------------
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# -----------------------------------------------------------------------------
# STEP 4: DADA2 denoising. This is long-running, so keep it inside a
# `screen` or `tmux` session on a remote server. Truncation lengths below
# (220/200) were chosen from this dataset's own quality plot; adjust based
# on your demux.qzv.
# -----------------------------------------------------------------------------
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 220 --p-trunc-len-r 200 \
  --p-n-threads 4 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoise-stats.qza \
  --verbose

# -----------------------------------------------------------------------------
# STEP 5: Summarize the feature table and denoising stats
# -----------------------------------------------------------------------------
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime metadata tabulate \
  --m-input-file denoise-stats.qza \
  --o-visualization denoise-stats.qzv

# -----------------------------------------------------------------------------
# STEP 6: Taxonomic classification against a SILVA classifier.
#
# Pre-trained classifiers must match your installed scikit-learn version.
# For qiime2-amplicon-2024.10 (scikit-learn 1.4.2), download the compatible
# full-length SILVA classifier below. It's about 500MB, so run this inside
# a screen or tmux session too:
#
#   wget -O silva-138-99-nb-classifier.qza \
#     "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
#
# A region-specific classifier (V4-only, for example) will usually be more
# accurate if one compatible with your QIIME2/scikit-learn version exists.
# The QIIME2 team stopped publishing pre-trained region-specific
# classifiers as of recent releases, so a full-length classifier is used
# here instead.
# -----------------------------------------------------------------------------
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-read-orientation same \
  --p-n-jobs 4 \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# -----------------------------------------------------------------------------
# STEP 7: Build the taxa barplot for exporting genus-level relative
# abundance tables (used downstream in the R analysis script)
# -----------------------------------------------------------------------------
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-barplot.qzv

echo "=== QIIME2 pipeline complete ==="
echo "Download demux.qzv, table.qzv, denoise-stats.qzv, taxonomy.qzv, and"
echo "taxa-barplot.qzv, then open them at https://view.qiime2.org to inspect"
echo "results and export the level-6 (genus) CSV for statistical analysis."
