#!/bin/bash
#SBATCH --job-name=phylum-by-code
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/barplot/all/phylum_by_code_all.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/barplot/all/phylum_by_code_all.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

WORK=~/sg/formal/merge4tubes
TABLE=$WORK/V4-merged-table-no-ch-filtered.qza
TAX=$WORK/V4-merged-taxonomy.qza
META=$WORK/merged_dolphin_4.txt
OUT=$WORK/barplot/all
mkdir -p "$OUT"

echo "[1/5] Collapse to Phylum & relative frequency…"
qiime taxa collapse \
  --i-table "$TABLE" \
  --i-taxonomy "$TAX" \
  --p-level 2 \
  --o-collapsed-table "$OUT/phylum-table.qza"

qiime feature-table relative-frequency \
  --i-table "$OUT/phylum-table.qza" \
  --o-relative-frequency-table "$OUT/phylum-relfreq.qza"

echo "[2/5] Export to TSV…"
rm -rf "$OUT/export"
qiime tools export --input-path "$OUT/phylum-relfreq.qza" --output-path "$OUT/export"
biom convert -i "$OUT/export/feature-table.biom" -o "$OUT/phylum_relab_raw.tsv" --to-tsv

echo "[3/5] Plot by code (with gaps)…"
python - <<'PY'
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser('~/sg/formal/merge4tubes')
OUT  = os.path.join(BASE, 'barplot', 'all')
meta_path = os.path.join(BASE, 'merged_dolphin_4.txt')
tsv_path  = os.path.join(OUT, 'phylum_relab_raw.tsv')

# metadata: detect columns
mraw = pd.read_csv(meta_path, sep='\t', dtype=str).fillna('')
sid_col = mraw.columns[0]
norm = {c: re.sub(r'[\s_]+','', c.strip().lower()) for c in mraw.columns}
code_col = next((k for k,v in norm.items() if v == 'code'), None)
if code_col is None: raise SystemExit("No 'code' column in metadata.")
meta = mraw[[sid_col, code_col]].rename(columns={sid_col:'SampleID', code_col:'code'})
meta['SampleID'] = meta['SampleID'].astype(str).str.strip()
meta['code'] = meta['code'].astype(str).str.strip()

# phylum table
df = pd.read_csv(tsv_path, sep='\t', skiprows=1)
df = df.rename(columns={'#OTU ID':'Phylum'})
def clean_phylum(x:str)->str:
    x = str(x)
    m = re.search(r'p__([^;]+)', x)
    return m.group(1) if m else x.replace('p__','').strip()
df['Phylum'] = df['Phylum'].map(clean_phylum)
df = df.set_index('Phylum')

# align & order samples by code then SampleID
sample_cols = [c for c in df.columns if c in set(meta['SampleID'])]
meta_ord = meta.loc[meta['SampleID'].isin(sample_cols)].sort_values(['code','SampleID'])
ordered_samples = meta_ord['SampleID'].tolist()
df = df[ordered_samples]

# top phyla + Other
TOP_N = 20
totals = df.sum(axis=1).sort_values(ascending=False)
top = totals.head(TOP_N).index.tolist()
top_df = df.loc[top]
other = (df.sum(axis=0) - top_df.sum(axis=0)).clip(lower=0)
plot_df = top_df.copy(); plot_df.loc['Other phyla'] = other
plot_df = plot_df.loc[list(totals.head(TOP_N).index) + ['Other phyla']]

# positions with gaps between code groups
gap = 0.5
xpos = []
current = 0.0
for code in sorted(meta_ord['code'].unique()):
    ns = (meta_ord['code']==code).sum()
    xpos.extend([current + i for i in range(ns)])
    current += ns + gap
xpos = np.array(xpos)

# colors
palette = list(plt.get_cmap('tab20').colors)
while len(palette) < plot_df.shape[0]:
    palette += list(plt.get_cmap('tab20b').colors)
colors = {phy: palette[i % len(palette)] for i, phy in enumerate(plot_df.index)}

# plot
fig, ax = plt.subplots(figsize=(18,6))
bottom = np.zeros(len(xpos))
for phy in plot_df.index[::-1]:
    vals = plot_df.loc[phy, :].values
    ax.bar(xpos, vals, width=0.9, bottom=bottom, color=colors[phy], linewidth=0)
    bottom += vals

# group separators & group labels
centers = []
start = 0
for code in sorted(meta_ord['code'].unique()):
    ns = (meta_ord['code']==code).sum()
    xs = xpos[start:start+ns]
    centers.append((xs.mean() if ns else start, code))
    start += ns + int(gap)
for i in range(len(centers)-1):
    ax.axvline((centers[i][0]+centers[i+1][0])/2, color='k', lw=0.5, alpha=0.4)

ax.set_xticks([c for c,_ in centers])
ax.set_xticklabels([c for _,c in centers], rotation=0, fontsize=10)

ax.set_ylabel('Relative abundance')
ax.set_ylim(0,1); ax.set_xlim(xpos.min()-0.6, xpos.max()+0.6)
ax.set_title('Bacterial phyla by sample — grouped by code', fontsize=14, pad=10)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.15)

# legend
handles = [plt.Rectangle((0,0),1,1,color=colors[p]) for p in plot_df.index]
ax.legend(handles, list(plot_df.index), bbox_to_anchor=(1.02,1), loc='upper left',
          frameon=False, fontsize=9)

# save
for ext in ('png','pdf','svg'):
    fig.savefig(os.path.join(OUT, f'phylum_by_code_barplot.{ext}'),
                dpi=300 if ext=='png' else None, bbox_inches='tight')
plt.close(fig)

# save order and table
plot_df.to_csv(os.path.join(OUT,'phylum_relab_top20_plus_other.tsv'), sep='\t')
meta_ord.to_csv(os.path.join(OUT,'sample_order_by_code.tsv'), sep='\t', index=False)
print("Done ->", OUT)
PY

echo "Finished."
