#!/bin/bash
#SBATCH --job-name=tube4o8-denoise-Oct       # Job name
#SBATCH --output=tube4o8-denoise-Oct.out     # Standard output log
#SBATCH --error=tube4o8-denoise-Oct.err      # Standard error log
#SBATCH --time=12:00:00                  # Time limit
#SBATCH --cpus-per-task=8                # Number of CPUs
#SBATCH --partition=normal               # Partition to use
#SBATCH --mail-type=BEGIN,END,FAIL       # Email notifications
#SBATCH --mail-user=shuai.gu@tamucc.edu  # Your email address

# Load QIIME 2 module
module load qiime2/2024.2


# Step 5: Denoise and merge reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/sg/formal/tube4/process/trimmed-demux-tube4o8.qza \
  --p-trunc-len-f 261 \
  --p-trunc-len-r 182 \
  --o-table ~/sg/formal/tube4/process/table-tube4o8.qza \
  --o-representative-sequences ~/sg/formal/tube4/process/rep-seqs-tube4o8.qza \
  --o-denoising-stats ~/sg/formal/tube4/process/denoising-stats-tube4o8.qza

# Step 6: Visualize outputs
qiime feature-table summarize \
  --i-table ~/sg/formal/tube4/process/table-tube4o8.qza \
  --o-visualization ~/sg/formal/tube4/process/table-tube4o8.qzv \
  --m-sample-metadata-file ~/sg/formal/tube4/tube4o8.txt

qiime feature-table tabulate-seqs \
  --i-data ~/sg/formal/tube4/process/rep-seqs-tube4o8.qza \
  --o-visualization ~/sg/formal/tube4/process/rep-seqs-tube4o8.qzv

qiime metadata tabulate \
  --m-input-file ~/sg/formal/tube4/process/denoising-stats-tube4o8.qza \
  --o-visualization ~/sg/formal/tube4/process/denoising-stats-tube4o8.qzv
