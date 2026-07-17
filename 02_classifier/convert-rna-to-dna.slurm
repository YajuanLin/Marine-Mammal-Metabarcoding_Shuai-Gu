#!/bin/bash
#SBATCH --job-name=convert-rna-to-dna
#SBATCH --output=convert-rna-to-dna.out
#SBATCH --error=convert-rna-to-dna.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

module load python

python convert_rna_to_dna.py
