#!/bin/bash
# =============================================================================
# 01_download_sra.sh
# Downloads raw 16S rRNA amplicon FASTQ files for Sun et al. 2025
# (BioProject PRJNA995545) from NCBI SRA.
#
# Requirements: sra-tools (prefetch, fasterq-dump), entrez-direct (esearch,
# efetch), gzip. All available via conda/bioconda.
#
# Usage:
#   chmod +x 01_download_sra.sh
#   ./01_download_sra.sh
# =============================================================================

set -e  # stop on first error

WORKDIR=~/sun2025_analysis
BIOPROJECT="PRJNA995545"

# -----------------------------------------------------------------------------
# Run this in a `screen` or `tmux` session first if working on a remote server
# -- the download loop below can take a long time and should survive an SSH
# disconnect:
#
#   screen -S sun_analysis
# -----------------------------------------------------------------------------

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
# Note: uses `for acc in $(cat ...)` rather than `while read acc; do ...`
# because prefetch/fasterq-dump inside a while-read loop can consume the
# loop's own stdin and cause it to exit early after only a few samples.
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

# Verify completeness: expect 140 files (70 samples x forward/reverse)
echo "File count:"
ls *.fastq.gz | wc -l

echo "Unique sample count:"
ls *.fastq.gz | sed 's/_[12]\.fastq\.gz//' | sort -u | wc -l
