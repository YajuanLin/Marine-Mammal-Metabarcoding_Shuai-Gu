#!/bin/bash
#SBATCH --job-name=alpha-rarefaction
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/alpha/alpha-rarefaction.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/alpha/alpha-rarefaction.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
ALPHA=$BASE/alpha
META=$BASE/merged_dolphin_4.txt
MAX_DEPTH=40000
STEPS=25

mkdir -p "$ALPHA" "$ALPHA/export"

echo "1) Rarefaction curves"
qiime diversity alpha-rarefaction \
  --i-table $BASE/V4-merged-table-no-ch-filtered.qza \
  --m-metadata-file $META \
  --p-max-depth $MAX_DEPTH \
  --p-steps $STEPS \
  --o-visualization $ALPHA/alpha-rarefaction_${MAX_DEPTH}.qzv

echo "2) Export sample counts"
qiime tools export \
  --input-path $BASE/V4-merged-table-no-ch-filtered.qza \
  --output-path $ALPHA/export

biom summarize-table \
  -i $ALPHA/export/feature-table.biom \
  -o $ALPHA/sample_counts.txt

echo "Done: $ALPHA/alpha-rarefaction_${MAX_DEPTH}.qzv  and  $ALPHA/sample_counts.txt"
