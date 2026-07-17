#!/bin/bash
#SBATCH --job-name=genus-collapse
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/genus_collapse.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/genus_collapse.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

# -----------------------------------------------------------------------------
# Genus-level (level 6) taxonomic collapse + relative frequency.
#
# This step was originally run interactively on 2025-10-08 (not saved as a
# script at the time). The command below is transcribed from the QIIME 2
# provenance embedded in genus-table.qza:
#   plugin: taxa   action: collapse   level: 6
#   inputs: table = V4-merged-table-no-ch-filtered.qza (uuid 35a670b0-...)
#           taxonomy = V4-merged-taxonomy.qza          (uuid 1922239a-...)
#   framework: QIIME 2 2024.2
# It is the genus-level sibling of barplot/all/phylum-by-code.slurm (level 2)
# and regenerates the genus-table.qza / genus-relfreq.qza consumed by the
# genus heatmap (Fig. 2) and the phylum+genus 2-panel (Fig. 5).
# -----------------------------------------------------------------------------

set -euo pipefail
module load qiime2/2024.2

WORK=~/sg/formal/merge4tubes
TABLE=$WORK/V4-merged-table-no-ch-filtered.qza
TAX=$WORK/V4-merged-taxonomy.qza
OUT=$WORK
mkdir -p "$OUT"

echo "[1/3] Collapse to Genus (level 6)…"
qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAX" \
  --p-level 6 \
  --o-collapsed-table "$OUT/genus-table.qza"

echo "[2/3] Relative frequency…"
qiime feature-table relative-frequency \
  --i-table "$OUT/genus-table.qza" \
  --o-relative-frequency-table "$OUT/genus-relfreq.qza"

echo "[3/3] Export to TSV…"
rm -rf "$OUT/genus-export"
qiime tools export --input-path "$OUT/genus-relfreq.qza" --output-path "$OUT/genus-export"
biom convert -i "$OUT/genus-export/feature-table.biom" -o "$OUT/genus_relab_raw.tsv" --to-tsv

echo "Finished."
