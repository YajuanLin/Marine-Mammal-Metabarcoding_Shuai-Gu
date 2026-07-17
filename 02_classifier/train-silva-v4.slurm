#!/bin/bash
#SBATCH --job-name=train-silva-v4
#SBATCH --output=train-silva-v4.out
#SBATCH --error=train-silva-v4.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load qiime2/2024.2

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-2-ref-seqs-V4.qza \
  --i-reference-taxonomy silva-138-2-tax.qza \
  --o-classifier ~/sg/ref/silva-138-2-classifier-V4.qza
