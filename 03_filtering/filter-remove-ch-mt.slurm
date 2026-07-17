#!/bin/bash
#SBATCH --job-name=filter-remove-ch-mt
#SBATCH --output=filter-remove-ch-mt.out
#SBATCH --error=filter-remove-ch-mt.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

WORK=~/sg/formal/merge4tubes
mkdir -p "$WORK"

# 1) Remove unassigned, chloroplast & mitochondria from rep-seqs
qiime taxa filter-seqs \
  --i-sequences       $WORK/V4-merged-rep-seqs.qza \
  --i-taxonomy        $WORK/V4-merged-taxonomy.qza \
  --p-exclude         Unassigned,Chloroplast,Mitochondria \
  --o-filtered-sequences $WORK/V4-merged-rep-seqs-no-ch.qza

# 2) Remove same taxa from feature table
qiime taxa filter-table \
  --i-table           $WORK/V4-merged-table.qza \
  --i-taxonomy        $WORK/V4-merged-taxonomy.qza \
  --p-exclude         Unassigned,Chloroplast,Mitochondria \
  --o-filtered-table  $WORK/V4-merged-table-no-ch.qza

# 3) Filter out any samples not in your metadata
qiime feature-table filter-samples \
  --i-table           $WORK/V4-merged-table-no-ch.qza \
  --m-metadata-file   $WORK/merged_dolphin_4.txt \
  --o-filtered-table  $WORK/V4-merged-table-no-ch-filtered.qza

# 3b) Sync representative sequences to the filtered table (keep only features present post-sample filtering)
qiime feature-table filter-seqs \
  --i-data               $WORK/V4-merged-rep-seqs-no-ch.qza \
  --i-table              $WORK/V4-merged-table-no-ch-filtered.qza \
  --o-filtered-data      $WORK/V4-merged-rep-seqs-no-ch-filtered.qza

# 4) Create the final barplot without unassigned/chloroplast/mito
qiime taxa barplot \
  --i-table           $WORK/V4-merged-table-no-ch-filtered.qza \
  --i-taxonomy        $WORK/V4-merged-taxonomy.qza \
  --m-metadata-file   $WORK/merged_dolphin_4.txt \
  --o-visualization   $WORK/taxa-barplot-V4-merged-no-ch-filtered.qzv
