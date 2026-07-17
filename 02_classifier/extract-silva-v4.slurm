#!/bin/bash
#SBATCH --job-name=extract-silva-v4
#SBATCH --output=extract-silva-v4.out
#SBATCH --error=extract-silva-v4.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load qiime2/2024.2

qiime feature-classifier extract-reads \
  --i-sequences silva-138-2-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 200 \
  --p-max-length 400 \
  --o-reads silva-138-2-ref-seqs-V4.qza
