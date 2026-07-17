#!/bin/bash
#SBATCH --job-name=import-silva-seqs
#SBATCH --output=import-silva-seqs.out
#SBATCH --error=import-silva-seqs.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load qiime2/2024.2

# Import reference sequences (uncompressed FASTA)
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path SILVA_138.2_SSURef_NR99_tax_silva_DNA.fasta \
  --output-path silva-138-2-seqs.qza
