#!/bin/bash
#SBATCH --job-name=tube2o8-denoise-July9th       # Job name
#SBATCH --output=tube2o8-denoise-July9th.out     # Standard output log
#SBATCH --error=tube2o8-denoise-July9th.err      # Standard error log
#SBATCH --time=12:00:00                  # Time limit
#SBATCH --cpus-per-task=8                # Number of CPUs
#SBATCH --partition=normal               # Partition to use
#SBATCH --mail-type=BEGIN,END,FAIL       # Email notifications
#SBATCH --mail-user=shuai.gu@tamucc.edu  # Your email address

# Load QIIME 2 module
module load qiime2/2024.2


# Step 5: Denoise and merge reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/sg/formal/tube2/process/trimmed-demux-tube2o8.qza \
  --p-trunc-len-f 272 \
  --p-trunc-len-r 184 \
  --o-table ~/sg/formal/tube2/process/table-tube2o8.qza \
  --o-representative-sequences ~/sg/formal/tube2/process/rep-seqs-tube2o8.qza \
  --o-denoising-stats ~/sg/formal/tube2/process/denoising-stats-tube2o8.qza

# Step 6: Visualize outputs
qiime feature-table summarize \
  --i-table ~/sg/formal/tube2/process/table-tube2o8.qza \
  --o-visualization ~/sg/formal/tube2/process/table-tube2o8.qzv \
  --m-sample-metadata-file ~/sg/formal/tube2/tube2o8.txt

qiime feature-table tabulate-seqs \
  --i-data ~/sg/formal/tube2/process/rep-seqs-tube2o8.qza \
  --o-visualization ~/sg/formal/tube2/process/rep-seqs-tube2o8.qzv

qiime metadata tabulate \
  --m-input-file ~/sg/formal/tube2/process/denoising-stats-tube2o8.qza \
  --o-visualization ~/sg/formal/tube2/process/denoising-stats-tube2o8.qzv
