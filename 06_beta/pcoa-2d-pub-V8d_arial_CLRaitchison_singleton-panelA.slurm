#!/bin/bash
#SBATCH --job-name=pcoa-panelA-CLRaitchison
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/beta/allbodysite/2D/%x.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/beta/allbodysite/2D/%x.err
#SBATCH --time=00:35:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
OUTDIR=$BASE/beta/allbodysite/2D
META=$BASE/merged_dolphin_4.txt
mkdir -p "$OUTDIR"

# ---- 1) Find table ----
TABLE_QZA=""
for CAND in \
  "$BASE/V4-merged-table-no-ch-filtered.qza" \
  "$BASE/beta/table-no-ch-filtered.qza" \
  "$BASE/beta/table.qza" \
  "$BASE/table.qza"
do
  if [ -f "$CAND" ]; then TABLE_QZA="$CAND"; break; fi
done

if [ -z "$TABLE_QZA" ]; then
  echo "ERROR: source table not found." >&2
  exit 1
fi

echo "Using source table:"
echo "  $TABLE_QZA"

# ---- 2) Apply SAME sample filters as original panel A ----
# Keep ONLY: genital_v, genital_p, mammary, gastric, rectal, blowhole
# Exclude: water samples and DolphinID 57
TAB_FILT=$OUTDIR/table_pub_counts_CLR_panelA.qza

qiime feature-table filter-samples \
  --i-table "$TABLE_QZA" \
  --m-metadata-file "$META" \
  --p-where "Type != 'water'
             AND DolphinID != '57'
             AND lower(BodySite) IN (
               'genital_v',
               'genital_p',
               'mammary',
               'gastric',
               'rectal',
               'blowhole'
             )" \
  --o-filtered-table "$TAB_FILT"

# ---- 3) Export filtered count table ----
CLRDIR=$OUTDIR/CLRaitchison_singletonsRemoved
EXPDIR_TABLE=$CLRDIR/export_table_counts
EXPDIR_ORD=$CLRDIR/export_ordination

rm -rf "$CLRDIR"
mkdir -p "$EXPDIR_TABLE" "$EXPDIR_ORD"

qiime tools export \
  --input-path "$TAB_FILT" \
  --output-path "$EXPDIR_TABLE"

biom convert \
  -i "$EXPDIR_TABLE/feature-table.biom" \
  -o "$EXPDIR_TABLE/feature-table.tsv" \
  --to-tsv

# ---- 4) Singleton removal + multiplicative replacement + CLR + Aitchison PCoA ----
python - <<'PY'
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa

# scikit-bio composition functions.
# QIIME2 versions may expose either multiplicative_replacement or multi_replace.
try:
    from skbio.stats.composition import clr, multiplicative_replacement
    zero_replace_func = multiplicative_replacement
    zero_replace_name = "multiplicative_replacement"
except ImportError:
    from skbio.stats.composition import clr, multi_replace
    zero_replace_func = multi_replace
    zero_replace_name = "multi_replace"

BASE = os.path.expanduser('~/sg/formal/merge4tubes')
OUTDIR = os.path.join(BASE, 'beta', 'allbodysite', '2D')
CLRDIR = os.path.join(OUTDIR, 'CLRaitchison_singletonsRemoved')
TABLE_TSV = os.path.join(CLRDIR, 'export_table_counts', 'feature-table.tsv')
EXPDIR_ORD = os.path.join(CLRDIR, 'export_ordination')
ORD = os.path.join(EXPDIR_ORD, 'ordination.txt')

os.makedirs(EXPDIR_ORD, exist_ok=True)

print("Using scikit-bio zero replacement function:", zero_replace_name)

# Robust BIOM TSV parse:
# line 1 is "# Constructed from biom file"
# line 2 is "#OTU ID <sample1> <sample2> ..."
tab = pd.read_csv(TABLE_TSV, sep='\t', skiprows=1, index_col=0)
tab.index.name = 'FeatureID'
tab.columns = tab.columns.astype(str)

# BIOM export is features x samples; transpose to samples x features.
counts = tab.T.apply(pd.to_numeric, errors='coerce').fillna(0)

# Basic checks
if counts.shape[0] == 0 or counts.shape[1] == 0:
    raise SystemExit("ERROR: count table has zero samples or zero features after parsing.")

if counts.index.nunique() != len(counts.index):
    dup = counts.index[counts.index.duplicated()].tolist()
    raise SystemExit(f"ERROR: duplicate sample IDs found: {dup[:10]}")

bad_ids = [s for s in counts.index if str(s).startswith('#')]
if bad_ids:
    raise SystemExit(f"ERROR: malformed sample IDs starting with '#': {bad_ids[:10]}")

# Remove all-zero features, just in case
counts = counts.loc[:, counts.sum(axis=0) > 0]

n_samples_before, n_features_before = counts.shape
zero_frac_before = (counts.values == 0).sum() / counts.size
sample_depth = counts.sum(axis=1)

# Remove singleton ASVs/features: retain features present in >=2 samples
feature_prevalence = (counts > 0).sum(axis=0)
min_prev = 2
keep_features = feature_prevalence >= min_prev
counts_filt = counts.loc[:, keep_features]

n_samples_after, n_features_after = counts_filt.shape
zero_frac_after = (counts_filt.values == 0).sum() / counts_filt.size

if n_features_after < 2:
    raise SystemExit("ERROR: fewer than 2 features retained after singleton removal.")

if (counts_filt.sum(axis=1) <= 0).any():
    bad_samples = counts_filt.index[counts_filt.sum(axis=1) <= 0].tolist()
    raise SystemExit(f"ERROR: samples with zero total counts after filtering: {bad_samples[:10]}")

# Convert to compositions / relative abundances before zero replacement
comp = counts_filt.div(counts_filt.sum(axis=1), axis=0)

# Multiplicative replacement for zeros
comp_repl = zero_replace_func(comp.values)
assert not np.any(comp_repl == 0), "zeros survived replacement; clr would produce -inf"

# CLR transform
clr_values = clr(comp_repl)
clr_df = pd.DataFrame(
    clr_values,
    index=counts_filt.index,
    columns=counts_filt.columns
)

# Euclidean distance on CLR-transformed composition = Aitchison distance
dist = squareform(pdist(clr_df.values, metric='euclidean'))
dm = DistanceMatrix(np.ascontiguousarray(dist, dtype=np.float64), ids=list(clr_df.index))

# PCoA
ord_res = pcoa(dm)
coords = ord_res.samples.copy()
eigvals = ord_res.eigvals.copy()
prop = ord_res.proportion_explained.copy()

# Save intermediate outputs
counts.to_csv(os.path.join(CLRDIR, 'counts_after_sample_filter.tsv'), sep='\t')
counts_filt.to_csv(os.path.join(CLRDIR, 'counts_after_singleton_removal.tsv'), sep='\t')
feature_prevalence.to_csv(os.path.join(CLRDIR, 'feature_prevalence_before_filter.tsv'), sep='\t', header=['prevalence'])
clr_df.to_csv(os.path.join(CLRDIR, 'clr_table.tsv'), sep='\t')
pd.DataFrame(dist, index=clr_df.index, columns=clr_df.index).to_csv(
    os.path.join(CLRDIR, 'aitchison_distance_matrix.tsv'),
    sep='\t'
)
coords.to_csv(os.path.join(EXPDIR_ORD, 'pcoa_coordinates.tsv'), sep='\t')
prop.to_csv(os.path.join(EXPDIR_ORD, 'proportion_explained.tsv'), sep='\t', header=False)

# Diagnostic summary
summary = {
    "n_samples": n_samples_before,
    "n_features_before_singleton_removal": n_features_before,
    "n_features_after_singleton_removal": n_features_after,
    "n_singleton_features_removed": int((~keep_features).sum()),
    "min_feature_prevalence_retained": min_prev,
    "sparsity_percent_before_singleton_removal": zero_frac_before * 100,
    "sparsity_percent_after_singleton_removal": zero_frac_after * 100,
    "min_sample_depth_after_sample_filter": sample_depth.min(),
    "median_sample_depth_after_sample_filter": sample_depth.median(),
    "mean_sample_depth_after_sample_filter": sample_depth.mean(),
    "max_sample_depth_after_sample_filter": sample_depth.max(),
    "zero_replacement_method": zero_replace_name,
    "transform": "CLR",
    "distance": "Euclidean on CLR-transformed compositions = Aitchison distance",
}
pd.Series(summary).to_csv(
    os.path.join(CLRDIR, 'CLR_Aitchison_summary.tsv'),
    sep='\t',
    header=False
)

# Write ordination.txt in a format compatible with your existing plotting parser
with open(ORD, 'w') as f:
    f.write("Eigvals\n")
    f.write("\t".join(str(v) for v in eigvals.values) + "\n\n")

    f.write("Proportion explained\n")
    f.write("\t".join(str(v) for v in prop.values) + "\n\n")

    f.write("Site\n")
    f.write(f"{coords.shape[0]} {coords.shape[1]}\n")
    for sid, row in coords.iterrows():
        vals = "\t".join(str(v) for v in row.values)
        f.write(f"{sid}\t{vals}\n")

print("CLR/Aitchison PCoA complete.")
print("--------------------------------")
for k, v in summary.items():
    print(f"{k}\t{v}")
print()
print("Wrote ordination:")
print(" ", ORD)
PY

# Clear matplotlib cache
rm -rf ~/.cache/matplotlib || true

# ---- 5) Plot panel A using same style as original, but corrected ordination ----
python - <<'PY'
import os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

BASE = os.path.expanduser('~/sg/formal/merge4tubes')
OUTDIR = os.path.join(BASE, 'beta', 'allbodysite', '2D')
META = os.path.join(BASE, 'merged_dolphin_4.txt')
CLRDIR = os.path.join(OUTDIR, 'CLRaitchison_singletonsRemoved')
ORD  = os.path.join(CLRDIR, 'export_ordination', 'ordination.txt')

# Enforce Arial
FONT_DIR = os.path.expanduser('~/sg/formal/merge4tubes/fonts/arial')
if os.path.isdir(FONT_DIR):
    for fn in os.listdir(FONT_DIR):
        p = os.path.join(FONT_DIR, fn)
        if p.lower().endswith('.ttf') and 'arial' in fn.lower():
            try:
                fm.fontManager.addfont(p)
            except Exception:
                pass
try:
    fm._rebuild()
except Exception:
    pass

mpl.rcParams.update({
    "font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","Liberation Sans","DejaVu Sans","Noto Sans","Nimbus Sans","sans-serif"],
    "pdf.fonttype":42,
    "ps.fonttype":42,
    "svg.fonttype":"none",
    "mathtext.fontset":"custom",
    "mathtext.rm":"Arial",
    "mathtext.it":"Arial:italic",
    "mathtext.bf":"Arial:bold",
    "axes.labelcolor":"black",
    "xtick.color":"black",
    "ytick.color":"black",
})

# Read metadata & map to canonical labels
import sys
sys.path.append('/home/sgu1/sg/formal/merge4tubes/palettes')
from palettes_pub import normalize_bodysite_series

meta = pd.read_csv(META, sep='\t', dtype=str).fillna('')
sid = meta.columns[0]
meta = meta.rename(columns={sid:'SampleID'})

meta['BodySite_pub'] = normalize_bodysite_series(
    meta.get('BodySite',''),
    meta.get('DolphinID','')
)

variant_map = {
    'genital_v':'Vaginal',
    'genital_p':'Penile',
    'mammary':'Mammary',
    'gastric':'Gastric',
    'rectal':'Rectal',
    'blowhole':'Respiratory',
}

bs_lower = meta.get('BodySite','').astype(str).str.lower()
mapped = bs_lower.map(variant_map)
meta['BodySite_pub'] = np.where(mapped.notna(), mapped, meta['BodySite_pub'])

BODY_ORDER = ['Vaginal','Penile','Mammary','Gastric','Rectal','Respiratory']
meta = meta[meta['BodySite_pub'].isin(BODY_ORDER)]

# Parse ordination.txt
text = open(ORD, 'r', encoding='utf-8', errors='ignore').read().splitlines()
float_re = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
pat_line = re.compile(r'^\S+\s+' + float_re + r'(?:\s+' + float_re + r')+')

def try_site_block(lines):
    rows = []
    i = 0
    hdr_re = re.compile(r'^(site|sites|sample|samples)\b', re.IGNORECASE)
    counts_re = re.compile(r'^\d+\s+\d+\s*$')
    while i < len(lines):
        s = lines[i].strip()
        if hdr_re.match(s):
            j = i + 1
            if j < len(lines) and counts_re.match(lines[j].strip()):
                j += 1
            while j < len(lines) and pat_line.match(lines[j].strip()):
                parts = lines[j].split()
                rows.append((parts[0], [float(x) for x in parts[1:]]))
                j += 1
            if rows:
                return rows
            i = j
        else:
            i += 1
    return []

def try_any_numeric_block(lines, valid):
    rows = []
    for s in lines:
        s = s.strip()
        if not pat_line.match(s):
            continue
        parts = s.split()
        sid = parts[0]
        if sid in valid:
            rows.append((sid, [float(x) for x in parts[1:]]))
    return rows

valid_ids = set(meta['SampleID'])
rows = try_site_block(text) or try_any_numeric_block(text, valid_ids)
if not rows:
    raise SystemExit("ERROR: No coordinates parsed from ordination.txt")

maxPC = max(len(v) for _, v in rows)
rec = []
for sid_, vec in rows:
    rec.append({
        'SampleID': sid_,
        **{f'PC{j+1}': (vec[j] if j < len(vec) else np.nan) for j in range(maxPC)}
    })

coords_df = pd.DataFrame(rec)
df = coords_df.merge(meta[['SampleID','BodySite_pub']], on='SampleID', how='inner')

# Hard check: do not silently drop samples during metadata merge
if df.shape[0] != coords_df.shape[0]:
    missing = sorted(set(coords_df['SampleID']) - set(df['SampleID']))
    raise SystemExit(
        "ERROR: Some ordination sample IDs did not merge with metadata. "
        f"Missing count: {len(missing)}. First few: {missing[:10]}"
    )

# Parse proportions
def parse_prop(path):
    try:
        L = [ln.strip() for ln in open(path, 'r', encoding='utf-8', errors='ignore')]
        for i, s in enumerate(L):
            if s.lower().startswith('proportion explained'):
                for x in L[i+1:]:
                    if x:
                        return [float(v) for v in x.split()]
    except Exception:
        pass
    return None

prop = parse_prop(ORD)

def axlabel(k):
    return f"PC{k} ({prop[k-1]*100:.2f}%)" if prop and len(prop) >= k else f"PC{k}"

# Same body-site colors as previous version
BODY_COLOR = {
    "Vaginal":     "#56B4E9",
    "Penile":      "#0072B2",
    "Mammary":     "#009E73",
    "Gastric":     "#E69F00",
    "Rectal":      "#CC79A7",
    "Respiratory": "#F0E442",
}

legend_handles = [
    Patch(facecolor=BODY_COLOR[bs], edgecolor='black', linewidth=0.6, label=bs)
    for bs in BODY_ORDER
]

def make_plot(a, b):
    x, y = f'PC{a}', f'PC{b}'
    if x not in df.columns or y not in df.columns:
        print(f"Skipping PC{a}-PC{b}: columns not found")
        return

    fig = plt.figure(figsize=(9.2, 5.2))
    ax = fig.add_subplot(111)

    for bs in BODY_ORDER:
        sub = df[df['BodySite_pub'] == bs]
        if sub.empty:
            continue
        ax.scatter(
            sub[x], sub[y],
            s=36,
            color=BODY_COLOR[bs],
            edgecolor='black',
            linewidths=0.4,
            alpha=0.95,
            label=bs
        )

    ax.set_xlabel(axlabel(a), fontsize=16)
    ax.set_ylabel(axlabel(b), fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add headroom above the cloud so the upper-left legend clears the points.
    # Cosmetic only (extends the view window; coordinates unchanged).
    ax.set_ylim(top=115)

    # Legend pinned to the upper-left corner (empty region of the PC1-PC2 panel)
    ax.legend(
        handles=legend_handles,
        title='BodySite',
        frameon=False,
        fontsize=15,
        title_fontsize=15,
        ncol=2,
        columnspacing=1.2,
        handletextpad=0.6,
        labelspacing=0.4,
        loc='upper left',
        bbox_to_anchor=(0.02, 0.98),
        bbox_transform=ax.transAxes,
    )

    stem = f"pcoa_axes{a}-{b}_BodySite_v8d_arial_CLRaitchison_singletonsRemoved_PC{a}-PC{b}"
    png = os.path.join(OUTDIR, f"{stem}.png")
    pdf = os.path.join(OUTDIR, f"{stem}.pdf")

    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')

    with open(os.path.join(OUTDIR, f"{stem}_font_report.txt"), "w", encoding="utf-8") as f:
        f.write("Panel A: CLR/Aitchison PCoA; singleton ASVs removed; multiplicative zero replacement; Arial; labels:16pt, ticks:14pt, legend:15pt; legend upper-left\n")

    plt.close(fig)
    print("Wrote:", png)
    print("Wrote:", pdf)

make_plot(1, 2)
make_plot(2, 3)
make_plot(1, 3)

print()
print("Merged ordination samples with metadata:", df.shape[0])
print("BodySite counts used for plotting:")
print(df['BodySite_pub'].value_counts().to_string())
PY

echo "Done."
echo "Corrected Panel A outputs are in:"
echo "  $OUTDIR"
echo
echo "Main PC1-PC2 file:"
echo "  $OUTDIR/pcoa_axes1-2_BodySite_v8d_arial_CLRaitchison_singletonsRemoved_PC1-PC2.png"
echo
echo "CLR/Aitchison diagnostic summary:"
echo "  $OUTDIR/CLRaitchison_singletonsRemoved/CLR_Aitchison_summary.tsv"
