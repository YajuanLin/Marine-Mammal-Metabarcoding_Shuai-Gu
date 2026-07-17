#!/bin/bash
#SBATCH --job-name=alpha-by-species
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-by-species.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-by-species.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
OUT=$BASE/alpha/more
DEPTH=10000
SPECIES_COL=Species
ID_COL=DolphinID
RAW_META=$BASE/merged_dolphin_4.txt
CLEAN_META=$OUT/merged_dolphin_4.clean.tsv

mkdir -p "$OUT"

# 0) Clean metadata: remove BOM/CR, trim header cell whitespace
if [[ ! -f "$RAW_META" ]]; then
  echo "ERROR: Metadata not found: $RAW_META"; exit 1
fi
# create cleaned copy
awk 'BEGIN{OFS="\t"} 
NR==1{
  gsub(/\r$/,""); 
  # split on tab; if your file is tab-delimited this preserves it
  n=split($0, a, "\t"); 
  for(i=1;i<=n;i++){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i]); }
  $0="";
  for(i=1;i<=n;i++){ $0 = $0 (i==1? "" : OFS) a[i]; }
  print; next
}
{ gsub(/\r$/,""); print }' "$RAW_META" > "$CLEAN_META"

# quick visibility
echo "Cleaned metadata saved to: $CLEAN_META"
echo "Header columns:"
head -n1 "$CLEAN_META" | tr '\t' '\n' | nl -ba

# robust header check (tab-delimited, trimmed)
hdr=$(head -n1 "$CLEAN_META")
for col in "$ID_COL" "$SPECIES_COL"; do
  if ! printf "%s" "$hdr" | tr '\t' '\n' | sed 's/^[[:space:]]\+//; s/[[:space:]]\+$//' | grep -Fxq "$col"; then
    echo "ERROR: Column '$col' not found after cleaning. Header is:"
    echo "$hdr"
    exit 1
  fi
done

echo "Validating metadata (writes $OUT/metadata.qzv)"
qiime metadata tabulate \
  --m-input-file "$CLEAN_META" \
  --o-visualization $OUT/metadata.qzv

echo "1) Remove ${ID_COL}=57"
qiime feature-table filter-samples \
  --i-table $BASE/V4-merged-table-no-ch-filtered.qza \
  --m-metadata-file "$CLEAN_META" \
  --p-where "${ID_COL} != '57'" \
  --o-filtered-table $OUT/table-no57.qza

for sp in Sb Tt; do
  echo "== species=${sp} =="
  SPDIR=$OUT/${sp}; mkdir -p "$SPDIR"

  echo "  A) Filter to ${SPECIES_COL}='${sp}'"
  qiime feature-table filter-samples \
    --i-table $OUT/table-no57.qza \
    --m-metadata-file "$CLEAN_META" \
    --p-where "${SPECIES_COL}='${sp}'" \
    --o-filtered-table $SPDIR/table-${sp}.qza

  echo "  B) Rarefy to ${DEPTH}"
  qiime feature-table rarefy \
    --i-table $SPDIR/table-${sp}.qza \
    --p-sampling-depth $DEPTH \
    --o-rarefied-table $SPDIR/rarefied-${DEPTH}.qza

  echo "  C) Alpha metrics"
  qiime diversity alpha --i-table $SPDIR/rarefied-${DEPTH}.qza --p-metric observed_features --o-alpha-diversity $SPDIR/observed_features.qza
  qiime diversity alpha --i-table $SPDIR/rarefied-${DEPTH}.qza --p-metric shannon            --o-alpha-diversity $SPDIR/shannon.qza
  qiime diversity alpha --i-table $SPDIR/rarefied-${DEPTH}.qza --p-metric simpson            --o-alpha-diversity $SPDIR/simpson.qza
  qiime diversity alpha --i-table $SPDIR/rarefied-${DEPTH}.qza --p-metric pielou_e           --o-alpha-diversity $SPDIR/pielou_e.qza
  qiime diversity alpha --i-table $SPDIR/rarefied-${DEPTH}.qza --p-metric chao1              --o-alpha-diversity $SPDIR/chao1.qza

  echo "  D) Export TSVs & merge"
  for m in observed_features shannon simpson pielou_e chao1; do
    expdir=$SPDIR/${m}-export
    rm -rf "$expdir"
    qiime tools export --input-path $SPDIR/${m}.qza --output-path "$expdir"
    mv "$expdir/alpha-diversity.tsv" $SPDIR/${m}.tsv
    rm -rf "$expdir"
  done

  awk -F'\t' '
    FNR==1 { next }
    FILENAME==ARGV[1] { ids[++n]=$1; obs[$1]=$2; next }
    FILENAME==ARGV[2] && FNR>1 { sha[$1]=$2; next }
    FILENAME==ARGV[3] && FNR>1 { sim[$1]=$2; next }
    FILENAME==ARGV[4] && FNR>1 { pie[$1]=$2; next }
    FILENAME==ARGV[5] && FNR>1 { cha[$1]=$2; next }
    END {
      print "SampleID\tobserved_features\tshannon\tsimpson\tpielou_e\tchao1";
      for (i=1;i<=n;i++){
        id=ids[i];
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", id, obs[id], sha[id], sim[id], pie[id], cha[id];
      }
    }
  ' $SPDIR/observed_features.tsv $SPDIR/shannon.tsv $SPDIR/simpson.tsv $SPDIR/pielou_e.tsv $SPDIR/chao1.tsv > $SPDIR/alpha_indices_${DEPTH}_${sp}.tsv

  echo "  -> $SPDIR/alpha_indices_${DEPTH}_${sp}.tsv"
done

echo "Done. Outputs in $OUT/Sb/ and $OUT/Tt/"
