#!/bin/bash
#SBATCH --job-name=tube1o8-July2nd       # Job name
#SBATCH --output=tube1o8-July2nd.out     # Standard output log
#SBATCH --error=tube1o8-July2nd.err      # Standard error log
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
  --input-path ~/sg/formal/tube1/process \
  --output-path ~/sg/formal/tube1/process/multiplexed-seqs-tube1.qza

# Step 2: Demultiplex
qiime cutadapt demux-paired \
  --i-seqs ~/sg/formal/tube1/process/multiplexed-seqs-tube1.qza \
  --m-forward-barcodes-file ~/sg/formal/tube1/tube1o8.txt \
  --m-forward-barcodes-column ForwardBarcode \
  --m-reverse-barcodes-file ~/sg/formal/tube1/tube1o8.txt \
  --m-reverse-barcodes-column ReverseBarcode \
  --o-per-sample-sequences ~/sg/formal/tube1/process/demux-tube1o8.qza \
  --o-untrimmed-sequences ~/sg/formal/tube1/process/untrimmed-tube1o8.qza

# Step 3: Summarize demultiplexed data
qiime demux summarize \
  --i-data ~/sg/formal/tube1/process/demux-tube1o8.qza \
  --o-visualization ~/sg/formal/tube1/process/demux-tube1o8.qzv

# Step 4: Trim primers/adapters
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ~/sg/formal/tube1/process/demux-tube1o8.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences ~/sg/formal/tube1/process/trimmed-demux-tube1o8.qza

qiime demux summarize \
  --i-data ~/sg/formal/tube1/process/trimmed-demux-tube1o8.qza \
  --o-visualization ~/sg/formal/tube1/process/trimmed-demux-tube1o8.qzv

# Step 5: Denoise and merge reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/sg/formal/tube1/process/trimmed-demux-tube1o8.qza \
  --p-trunc-len-f 244 \
  --p-trunc-len-r 198 \
  --o-table ~/sg/formal/tube1/process/table-tube1o8.qza \
  --o-representative-sequences ~/sg/formal/tube1/process/rep-seqs-tube1o8.qza \
  --o-denoising-stats ~/sg/formal/tube1/process/denoising-stats-tube1o8.qza

# Step 6: Visualize outputs
qiime feature-table summarize \
  --i-table ~/sg/formal/tube1/process/table-tube1o8.qza \
  --o-visualization ~/sg/formal/tube1/process/table-tube1o8.qzv \
  --m-sample-metadata-file ~/sg/formal/tube1/tube1o8.txt

qiime feature-table tabulate-seqs \
  --i-data ~/sg/formal/tube1/process/rep-seqs-tube1o8.qza \
  --o-visualization ~/sg/formal/tube1/process/rep-seqs-tube1o8.qzv

qiime metadata tabulate \
  --m-input-file ~/sg/formal/tube1/process/denoising-stats-tube1o8.qza \
  --o-visualization ~/sg/formal/tube1/process/denoising-stats-tube1o8.qzv
