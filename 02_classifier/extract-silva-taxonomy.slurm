#!/bin/bash
#SBATCH --job-name=extract-silva-taxonomy
#SBATCH --output=extract-silva-taxonomy.out
#SBATCH --error=extract-silva-taxonomy.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load python
pip install --user biopython

python extract_silva_taxonomy.py
