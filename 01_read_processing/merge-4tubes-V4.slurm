#!/bin/bash
#SBATCH --job-name=merge-4tubes-V4
#SBATCH --output=merge-4tubes-V4.out
#SBATCH --error=merge-4tubes-V4.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load qiime2/2024.2

WORK=~/sg/formal/merge4tubes
mkdir -p "$WORK"

# 1) Merge the four feature tables (add tube4)
qiime feature-table merge \
  --i-tables ~/sg/formal/tube1/process/table-tube1o8.qza \
  --i-tables ~/sg/formal/tube2/process/table-tube2o8.qza \
  --i-tables ~/sg/formal/tube3/process/table-tube3o8.qza \
  --i-tables ~/sg/formal/tube4/process/table-tube4o8.qza \
  --o-merged-table $WORK/V4-merged-table.qza

# 2) Merge rep‐seqs (add tube4) and classify merged sequences
qiime feature-table merge-seqs \
  --i-data ~/sg/formal/tube1/process/rep-seqs-tube1o8.qza \
  --i-data ~/sg/formal/tube2/process/rep-seqs-tube2o8.qza \
  --i-data ~/sg/formal/tube3/process/rep-seqs-tube3o8.qza \
  --i-data ~/sg/formal/tube4/process/rep-seqs-tube4o8.qza \
  --o-merged-data $WORK/V4-merged-rep-seqs.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ~/sg/ref/silva-138-2-classifier-V4.qza \
  --i-reads $WORK/V4-merged-rep-seqs.qza \
  --o-classification $WORK/V4-merged-taxonomy.qza \
  --p-n-jobs 8
