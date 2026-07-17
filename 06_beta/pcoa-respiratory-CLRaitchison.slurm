#!/bin/bash
#SBATCH --job-name=pcoa-respiratory-CLRaitchison
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/beta/respiratory/2D/%x.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/beta/respiratory/2D/%x.err
#SBATCH --time=00:35:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
GEN=$BASE/beta/respiratory/2D
META=$BASE/merged_dolphin_4.txt
mkdir -p "$GEN"

# Per-site config
SITE_NAME="respiratory"
CLRDIR=$GEN/CLRaitchison_singletonsRemoved
EXPDIR_TABLE=$CLRDIR/export_table_counts
EXPDIR_ORD=$CLRDIR/export_ordination

# ---- 1) Find source table ----
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
echo "Using source table: $TABLE_QZA"

# ---- 2) Filter samples to this site ----
# respiratory. Same exclusions as all panels.
TAB_FILT=$GEN/table_${SITE_NAME}_counts_CLR.qza
qiime feature-table filter-samples \
  --i-table "$TABLE_QZA" \
  --m-metadata-file "$META" \
  --p-where "Type != 'water'
             AND DolphinID != '57'
             AND lower(BodySite) IN ('blowhole')" \
  --o-filtered-table "$TAB_FILT"

# ---- 3) Export filtered count table ----
rm -rf "$CLRDIR"
mkdir -p "$EXPDIR_TABLE" "$EXPDIR_ORD"

qiime tools export \
  --input-path "$TAB_FILT" \
  --output-path "$EXPDIR_TABLE"

biom convert \
  -i "$EXPDIR_TABLE/feature-table.biom" \
  -o "$EXPDIR_TABLE/feature-table.tsv" \
  --to-tsv

# ---- 4) Prevalence filter + multiplicative replacement + CLR + Aitchison PCoA ----
export GEN CLRDIR SITE_NAME
python - <<'PY'
import os
import math
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa

try:
    from skbio.stats.composition import clr, multiplicative_replacement
    zero_replace_func = multiplicative_replacement
    zero_replace_name = "multiplicative_replacement"
except ImportError:
    from skbio.stats.composition import clr, multi_replace
    zero_replace_func = multi_replace
    zero_replace_name = "multi_replace"

CLRDIR = os.environ['CLRDIR']
SITE_NAME = os.environ['SITE_NAME']
TABLE_TSV = os.path.join(CLRDIR, 'export_table_counts', 'feature-table.tsv')
EXPDIR_ORD = os.path.join(CLRDIR, 'export_ordination')
ORD = os.path.join(EXPDIR_ORD, 'ordination.txt')
os.makedirs(EXPDIR_ORD, exist_ok=True)

print(f"[{SITE_NAME}] zero replacement function:", zero_replace_name)

# Robust BIOM TSV parse
tab = pd.read_csv(TABLE_TSV, sep='\t', skiprows=1, index_col=0)
tab.index.name = 'FeatureID'
tab.columns = tab.columns.astype(str)
counts = tab.T.apply(pd.to_numeric, errors='coerce').fillna(0)

if counts.shape[0] == 0 or counts.shape[1] == 0:
    raise SystemExit("ERROR: empty count table after parsing.")
if counts.index.nunique() != len(counts.index):
    dup = counts.index[counts.index.duplicated()].tolist()
    raise SystemExit(f"ERROR: duplicate sample IDs: {dup[:10]}")
bad_ids = [s for s in counts.index if str(s).startswith('#')]
if bad_ids:
    raise SystemExit(f"ERROR: malformed sample IDs starting with '#': {bad_ids[:10]}")

# Drop all-zero features
counts = counts.loc[:, counts.sum(axis=0) > 0]
n_samples_before, n_features_before = counts.shape
zero_frac_before = (counts.values == 0).sum() / counts.size
sample_depth = counts.sum(axis=1)

# Remove singleton ASVs/features: retain features present in at least 2 samples.
min_prev = 2
feature_prevalence = (counts > 0).sum(axis=0)
keep_features = feature_prevalence >= min_prev
counts_filt = counts.loc[:, keep_features]

n_samples_after, n_features_after = counts_filt.shape
zero_frac_after = (counts_filt.values == 0).sum() / counts_filt.size

if n_features_after < 2:
    raise SystemExit("ERROR: <2 features retained after prevalence filter.")
if (counts_filt.sum(axis=1) <= 0).any():
    bad = counts_filt.index[counts_filt.sum(axis=1) <= 0].tolist()
    raise SystemExit(f"ERROR: samples with zero counts after filtering: {bad[:10]}")

# Closure -> multiplicative replacement -> CLR
comp = counts_filt.div(counts_filt.sum(axis=1), axis=0)
comp_repl = zero_replace_func(comp.values)
assert not np.any(comp_repl == 0), "zeros survived replacement; clr would produce -inf"
clr_values = clr(comp_repl)
clr_df = pd.DataFrame(clr_values, index=counts_filt.index, columns=counts_filt.columns)

# Aitchison distance = Euclidean on CLR
dist = squareform(pdist(clr_df.values, metric='euclidean'))
dm = DistanceMatrix(np.ascontiguousarray(dist, dtype=np.float64), ids=list(clr_df.index))

# PCoA
ord_res = pcoa(dm)
coords = ord_res.samples.copy()
eigvals = ord_res.eigvals.copy()
prop = ord_res.proportion_explained.copy()

# Save intermediates
clr_df.to_csv(os.path.join(CLRDIR, 'clr_table.tsv'), sep='\t')
pd.DataFrame(dist, index=clr_df.index, columns=clr_df.index).to_csv(
    os.path.join(CLRDIR, 'aitchison_distance_matrix.tsv'), sep='\t')
coords.to_csv(os.path.join(EXPDIR_ORD, 'pcoa_coordinates.tsv'), sep='\t')
prop.to_csv(os.path.join(EXPDIR_ORD, 'proportion_explained.tsv'), sep='\t', header=False)

summary = {
    "site": SITE_NAME,
    "n_samples": n_samples_before,
    "n_features_before_prevalence_filter": n_features_before,
    "n_features_after_prevalence_filter": n_features_after,
    "n_features_removed": int((~keep_features).sum()),
    "min_feature_prevalence_retained": min_prev,
    "prevalence_rule": "singletons removed (min prevalence = 2 samples)",
    "sparsity_percent_before": zero_frac_before * 100,
    "sparsity_percent_after": zero_frac_after * 100,
    "min_sample_depth": sample_depth.min(),
    "median_sample_depth": sample_depth.median(),
    "max_sample_depth": sample_depth.max(),
    "zero_replacement_method": zero_replace_name,
    "transform": "CLR",
    "distance": "Euclidean on CLR (Aitchison)",
}
pd.Series(summary).to_csv(os.path.join(CLRDIR, 'CLR_Aitchison_summary.tsv'), sep='\t', header=False)

with open(ORD, 'w') as f:
    f.write("Eigvals\n")
    f.write("\t".join(str(v) for v in eigvals.values) + "\n\n")
    f.write("Proportion explained\n")
    f.write("\t".join(str(v) for v in prop.values) + "\n\n")
    f.write("Site\n")
    f.write(f"{coords.shape[0]} {coords.shape[1]}\n")
    for sid, row in coords.iterrows():
        f.write(f"{sid}\t" + "\t".join(str(v) for v in row.values) + "\n")

print(f"[{SITE_NAME}] CLR/Aitchison PCoA complete.")
for k, v in summary.items():
    print(f"  {k}\t{v}")
print(f"[{SITE_NAME}] wrote ordination: {ORD}")
PY

# Clear matplotlib cache
rm -rf ~/.cache/matplotlib || true

# ---- 5) Plot: facility color + species/sex marker, legend-less (same scheme as old per-panel) ----
export GEN BASE CLRDIR SITE_NAME
python - <<'PY'
import os, re, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib as mpl
from matplotlib import font_manager as fm
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.text import Text
from collections import Counter

BASE = os.environ['BASE']
GEN  = os.environ['GEN']
CLRDIR = os.environ['CLRDIR']
SITE_NAME = os.environ['SITE_NAME']
VERSION = "v8d_arial_CLRaitchison_singletonsRemoved"
ord_path  = os.path.join(CLRDIR, 'export_ordination', 'ordination.txt')
meta_path = os.path.join(BASE, 'merged_dolphin_4.txt')
coord_tsv = os.path.join(GEN, f'pcoa_coords_{SITE_NAME}_gender_code_{VERSION}.tsv')

# ====== Force ARIAL ======
FONT_DIR = os.path.expanduser("~/sg/formal/merge4tubes/fonts/arial")
FOUND = []
if os.path.isdir(FONT_DIR):
    for fn in os.listdir(FONT_DIR):
        p = os.path.join(FONT_DIR, fn)
        if p.lower().endswith(".ttf") and "arial" in fn.lower():
            try:
                fm.fontManager.addfont(p); FOUND.append(p)
            except Exception:
                pass
try: fm._rebuild()
except Exception: pass

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial","Helvetica","Liberation Sans","DejaVu Sans","Noto Sans","Nimbus Sans","sans-serif"],
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

def export_font_report(fig, path, extra=None):
    lines=[f"Matplotlib version: {matplotlib.__version__}",
           f"rcParams[font.family]: {mpl.rcParams.get('font.family')}",
           f"User-registered Arial files: {len(FOUND)}"]
    resolved = findfont(FontProperties(family=['Arial','sans-serif']))
    lines.append(f"Resolved 'Arial' path: {resolved}\n")
    names=[]
    for i,t in enumerate(fig.findobj(match=Text),1):
        fp=t.get_fontproperties(); names.append(fp.get_name())
    for n,c in Counter(names).most_common():
        lines.append(f"{n}: {c}")
    if extra:
        for k,v in extra.items(): lines.append(f"{k}: {v}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

sys.path.append('/home/sgu1/sg/formal/merge4tubes/palettes')
from palettes_pub import normalize_bodysite_series

# Facility colors (Paul Tol Bright) -- identical to old per-panel scripts
CODE_ORDER = ["HI1", "HI2", "FL1", "FL2", "FL3", "FL4"]
CODE_COLOR = {
    "HI1": "#4477AA", "HI2": "#EE6677", "FL1": "#228833",
    "FL2": "#AA3377", "FL3": "#66CCEE", "FL4": "#CCBB44",
    "Unknown": "#999999",
}

def load_meta(path):
    m = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    def pick(cols, *alts):
        low = {c.lower().replace(' ','').replace('_',''): c for c in cols}
        for a in alts:
            k = a.lower().replace(' ','').replace('_','')
            if k in low: return low[k]
        return None
    sid = m.columns[0]
    code = pick(m.columns, 'code')
    bs   = pick(m.columns, 'bodysite','body_site','bodysites')
    sp   = pick(m.columns, 'species')
    did  = pick(m.columns, 'dolphinid','animalid','id')
    typ  = pick(m.columns, 'type','sampletype')
    sex  = pick(m.columns, 'gender','sex')
    m = m.rename(columns={sid:'SampleID', code:'Code', bs:'BodySite',
                          sp:'Species', did:'DolphinID', typ:'Type', sex:'Gender'})
    for c in ['SampleID','Code','BodySite','Species','DolphinID','Type','Gender']:
        if c in m.columns:
            m[c] = m[c].astype(str).str.strip()
        else:
            m[c] = ''
    m['BodySite_pub'] = normalize_bodysite_series(m['BodySite'], m['DolphinID'])
    g = m['Gender'].str.lower()
    m['is_female'] = g.isin({'f','female'})
    m['is_male']   = g.isin({'m','male'})
    return m

def parse_prop(path):
    try:
        lines=[ln.strip() for ln in open(path,'r',encoding='utf-8',errors='ignore')]
        for i,s in enumerate(lines):
            if s.lower().startswith('proportion explained'):
                for x in lines[i+1:]:
                    if x:
                        return [float(v) for v in x.split()]
    except Exception:
        pass
    return None

def axis_label(k, prop):
    return f"PC{k} ({prop[k-1]*100:.2f}%)" if (prop and len(prop)>=k) else f"PC{k}"

# ---- metadata filters: SAME as old genital panel (effective set Vaginal+Penile) ----
meta = load_meta(meta_path)
meta = meta[meta['DolphinID']!='57']
meta = meta[meta['Type'].str.lower()!='water']
meta = meta[~meta['BodySite_pub'].isin({'Skin','Milk','Maternal wall','Fetal wall'})]
meta = meta[meta['BodySite_pub']=='Respiratory']
if meta.empty:
    raise SystemExit(f"No {SITE_NAME} samples after filtering.")
valid_ids = set(meta['SampleID'])

# ---- Parse ordination.txt ----
text = open(ord_path, 'r', encoding='utf-8', errors='ignore').read().splitlines()
float_re = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
pat_line = re.compile(r'^\S+\s+' + float_re + r'(?:\s+' + float_re + r')+')

def try_site_block(lines):
    rows=[]; i=0
    hdr_re = re.compile(r'^(site|sites|sample|samples)\b', re.IGNORECASE)
    counts_re = re.compile(r'^\d+\s+\d+\s*$')
    while i < len(lines):
        s = lines[i].strip()
        if hdr_re.match(s):
            j = i+1
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
    rows=[]
    for s in lines:
        s = s.strip()
        if not pat_line.match(s): continue
        parts = s.split()
        if parts[0] in valid:
            rows.append((parts[0], [float(x) for x in parts[1:]]))
    return rows

rows = try_site_block(text) or try_any_numeric_block(text, valid_ids)
if not rows:
    raise SystemExit("Could not locate coordinates in ordination.txt")

maxPC = max(len(v) for _, v in rows)
rec = [{'SampleID': sid, **{f'PC{j+1}': (vec[j] if j < len(vec) else np.nan)
        for j in range(maxPC)}} for sid, vec in rows]
coords_df = pd.DataFrame(rec)
df = coords_df.merge(
    meta[['SampleID','BodySite_pub','Species','Gender','Code','is_female','is_male']],
    on='SampleID', how='inner')

# Hard check: no silent sample loss in merge
if df.shape[0] != coords_df.shape[0]:
    missing = sorted(set(coords_df['SampleID']) - set(df['SampleID']))
    raise SystemExit(f"ERROR: {len(missing)} ordination IDs missing from metadata. "
                     f"First few: {missing[:10]}")

for pc in [c for c in ['PC1','PC2','PC3'] if c in df.columns]:
    df[pc] = pd.to_numeric(df[pc], errors='coerce')

keep=['SampleID','BodySite_pub','Species','Gender','Code','is_female','is_male','PC1','PC2','PC3']
df[[c for c in keep if c in df.columns]].to_csv(coord_tsv, sep='\t', index=False)

prop = parse_prop(ord_path)
present_codes = [c for c in CODE_ORDER if c in set(df['Code'])]
code_color = {c: CODE_COLOR.get(c,'#999999') for c in present_codes}

is_Tt = df['Species'].str.upper().str.startswith('TT')
is_Sb = df['Species'].str.upper().str.startswith('SB')
female = df['is_female'].astype(bool)
male   = df['is_male'].astype(bool)
tt_f = df[ is_Tt & female].copy()   # solid circle
tt_m = df[ is_Tt & male].copy()     # hollow circle (code-colored edge)
sb_f = df[ is_Sb & female].copy()   # solid square

def scatter(ax, dsub, x, y, marker='o', filled=True, outline_by_code=False):
    if dsub.empty: return
    cols = dsub['Code'].map(lambda c: code_color.get(c,'#999999'))
    if filled:
        ax.scatter(dsub[x], dsub[y], c=cols, s=36, alpha=0.95, marker=marker, edgecolor='none')
    else:
        edgecols = cols if outline_by_code else 'k'
        ax.scatter(dsub[x], dsub[y], facecolors='none', s=36, alpha=0.95,
                   marker=marker, edgecolor=edgecols, linewidths=1.0)

def make_plot(a, b):
    x, y = f'PC{a}', f'PC{b}'
    if x not in df.columns or y not in df.columns:
        return
    fig = plt.figure(figsize=(9.2,5.2))
    ax  = fig.add_subplot(111)
    scatter(ax, tt_f, x, y, marker='o', filled=True)
    scatter(ax, tt_m, x, y, marker='o', filled=False, outline_by_code=True)
    scatter(ax, sb_f, x, y, marker='s', filled=True)
    ax.set_xlabel(axis_label(a, prop), fontsize=16)
    ax.set_ylabel(axis_label(b, prop), fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(False)
    stem = f"pcoa_axes{a}-{b}_{SITE_NAME}_genderFacility_{VERSION}_PC{a}-PC{b}"
    png = os.path.join(GEN, f"{stem}.png")
    pdf = os.path.join(GEN, f"{stem}.pdf")
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    export_font_report(fig, os.path.join(GEN, f"{stem}_font_report.txt"),
                       {"output_files": f"{os.path.basename(png)}, {os.path.basename(pdf)}"})
    plt.close(fig)
    print("Wrote:", png)

make_plot(1,2)
make_plot(2,3)
make_plot(1,3)

print()
print(f"[{SITE_NAME}] merged samples plotted:", df.shape[0])
print(f"[{SITE_NAME}] facility (Code) counts:")
print(df['Code'].value_counts().to_string())
PY

echo "Done. Panel ${SITE_NAME} CLR/Aitchison outputs in: $GEN"
echo "Main PC1-PC2:"
echo "  $GEN/pcoa_axes1-2_${SITE_NAME}_genderFacility_v8d_arial_CLRaitchison_singletonsRemoved_PC1-PC2.png"
