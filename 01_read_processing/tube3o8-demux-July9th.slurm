#!/bin/bash
#SBATCH --job-name=tube3o8-July9th       # Job name
#SBATCH --output=tube3o8-July9th.out     # Standard output log
#SBATCH --error=tube3o8-July9th.err      # Standard error log
#SBATCH --time=12:00:00                  # Time limit
#SBATCH --cpus-per-task=8                # Number of CPUs
#SBATCH --partition=normal               # Partition to use
#SBATCH --mail-type=BEGIN,END,FAIL       # Email notifications
#SBATCH --mail-user=shuai.gu@tamucc.edu  # Your email address

# Load QIIME 2 module
module load qiime2/2024.2

# Step 1: Import
qiime tools import \
  --type 'MultiplexedPairedEndBarcodeInSequence' \
  --input-path ~/sg/formal/tube3/process \
  --output-path ~/sg/formal/tube3/process/multiplexed-seqs-tube3.qza

# Step 2: Demultiplex
qiime cutadapt demux-paired \
  --i-seqs ~/sg/formal/tube3/process/multiplexed-seqs-tube3.qza \
  --m-forward-barcodes-file ~/sg/formal/tube3/tube3o8.txt \
  --m-forward-barcodes-column ForwardBarcode \
  --m-reverse-barcodes-file ~/sg/formal/tube3/tube3o8.txt \
  --m-reverse-barcodes-column ReverseBarcode \
  --o-per-sample-sequences ~/sg/formal/tube3/process/demux-tube3o8.qza \
  --o-untrimmed-sequences ~/sg/formal/tube3/process/untrimmed-tube3o8.qza

# Step 3: Summarize demultiplexed data
qiime demux summarize \
  --i-data ~/sg/formal/tube3/process/demux-tube3o8.qza \
  --o-visualization ~/sg/formal/tube3/process/demux-tube3o8.qzv

# Step 4: Trim primers/adapters
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ~/sg/formal/tube3/process/demux-tube3o8.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences ~/sg/formal/tube3/process/trimmed-demux-tube3o8.qza

qiime demux summarize \
  --i-data ~/sg/formal/tube3/process/trimmed-demux-tube3o8.qza \
  --o-visualization ~/sg/formal/tube3/process/trimmed-demux-tube3o8.qzv

