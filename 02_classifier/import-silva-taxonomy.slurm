#!/bin/bash
#SBATCH --job-name=import-silva-taxonomy
#SBATCH --output=import-silva-taxonomy.out
#SBATCH --error=import-silva-taxonomy.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load qiime2/2024.2

# Import SILVA taxonomy (non-prefixed ranks are acceptable)
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path silva-138-2-taxonomy.txt \
  --output-path silva-138-2-tax.qza
