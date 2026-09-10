#!/bin/bash
# =============================================================================
# 01_download_sra.sh
# Downloads raw 16S rRNA amplicon FASTQ files for Sun et al. 2025
# (BioProject PRJNA995545) from NCBI SRA.
#
# Requirements: sra-tools (prefetch, fasterq-dump), entrez-direct (esearch,
# efetch), gzip. All available via conda/bioconda.
#
# -----------------------------------------------------------------------------
# Connecting to a remote server (skip this if running locally)
#
# This script was originally run on a remote Linux server over SSH. If
# you're doing the same, connect first and activate your QIIME2 environment
# before running anything below:
#
#   ssh -p <port> <username>@<server-address>
#   conda activate qiime2-amplicon-2024.10
#
# Long-running steps in this script (the download loop) should survive an
# SSH disconnect, so start a screen session before continuing:
#
#   screen -S sun_analysis
#
# If you get disconnected and reconnect later, resume the same session
# instead of starting a new one:
#
#   screen -r sun_analysis
# -----------------------------------------------------------------------------
#
# Usage:
#   chmod +x 01_download_sra.sh
#   ./01_download_sra.sh
# =============================================================================

set -e  # stop on first error

WORKDIR=~/sun2025_analysis
BIOPROJECT="PRJNA995545"

mkdir -p "$WORKDIR" && cd "$WORKDIR"

# Fetch the run list for the BioProject and extract accession IDs
esearch -db sra -query "$BIOPROJECT" | efetch -format runinfo > runinfo.csv
cut -d',' -f1 runinfo.csv | tail -n +2 > accession_list.txt

echo "Sample count:"
wc -l accession_list.txt
# Expected: 70 (matches the paper's 70 skin samples)

# -----------------------------------------------------------------------------
# Download and convert each sample to gzipped paired-end FASTQ.
#
# Note: this uses `for acc in $(cat ...)` rather than `while read acc; do
# ...`. A while-read loop can have its stdin consumed by prefetch or
# fasterq-dump running inside it, which makes the loop exit early after
# only a few samples.
# -----------------------------------------------------------------------------
mkdir -p raw_fastq && cd raw_fastq

for acc in $(cat ../accession_list.txt); do
  echo "=== Downloading $acc ==="
  prefetch "$acc" -O . --max-size 1g
  fasterq-dump "$acc" --split-files -O . -e 4
  gzip -f "${acc}"*.fastq
  rm -rf "$acc"
done

echo "=== ALL DOWNLOADS COMPLETE ==="

# Verify completeness. Expect 140 files: 70 samples, forward and reverse each.
echo "File count:"
ls *.fastq.gz | wc -l

echo "Unique sample count:"
ls *.fastq.gz | sed 's/_[12]\.fastq\.gz//' | sort -u | wc -l
