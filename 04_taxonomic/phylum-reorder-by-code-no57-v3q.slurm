#!/bin/bash
#SBATCH --job-name=phylum-reorder-by-code-no57-v3q
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/barplot/all/phylum_reorder_by_code_no57_v3q.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/barplot/all/phylum_reorder_by_code_no57_v3q.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
OUT=$BASE/barplot/all
META=$BASE/merged_dolphin_4.txt
TSV=$OUT/phylum_relab_raw.tsv
RELQZA=$OUT/phylum-relfreq.qza
PHYLUMQZA=$OUT/phylum-table.qza

# Ensure rel-abundance TSV exists
if [ ! -f "$TSV" ]; then
  if [ -f "$RELQZA" ]; then
    rm -rf "$OUT/export"
    qiime tools export --input-path "$RELQZA" --output-path "$OUT/export"
    biom convert -i "$OUT/export/feature-table.biom" -o "$TSV" --to-tsv
  elif [ -f "$PHYLUMQZA" ]; then
    qiime feature-table relative-frequency --i-table "$PHYLUMQZA" --o-relative-frequency-table "$RELQZA"
    rm -rf "$OUT/export"
    qiime tools export --input-path "$RELQZA" --output-path "$OUT/export"
    biom convert -i "$OUT/export/feature-table.biom" -o "$TSV" --to-tsv
  else
    echo "Missing phylum inputs (phylum-table/relfreq)." >&2; exit 1
  fi
fi

# Clear per-user Matplotlib cache to avoid stale font fallback
rm -rf ~/.cache/matplotlib || true

python - <<'PY'
import os, re, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties, findfont
from collections import Counter

# ---------- Register user-provided Arial fonts ----------
FONT_DIR = os.path.expanduser("~/sg/formal/merge4tubes/fonts/arial")
FOUND = []
if os.path.isdir(FONT_DIR):
    for fn in os.listdir(FONT_DIR):
        p = os.path.join(FONT_DIR, fn)
        if p.lower().endswith(".ttf") and "arial" in fn.lower():
            try:
                font_manager.fontManager.addfont(p)
                FOUND.append(p)
            except Exception:
                pass
try:
    font_manager._rebuild()
except Exception:
    pass

# ---------- Force Arial across all text ----------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial","Helvetica","Liberation Sans","DejaVu Sans","Noto Sans","Nimbus Sans","sans-serif"],
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})

sys.path.append('/home/sgu1/sg/formal/merge4tubes/palettes')
from palettes_pub import normalize_bodysite_series

BASE = os.path.expanduser('~/sg/formal/merge4tubes')
OUT  = os.path.join(BASE, 'barplot', 'all')
meta_path = os.path.join(BASE, 'merged_dolphin_4.txt')
tsv_path  = os.path.join(OUT, 'phylum_relab_raw.tsv')

# --- metadata ---
mraw = pd.read_csv(meta_path, sep='\t', dtype=str).fillna('')
sid_col = mraw.columns[0]
norm = {c: re.sub(r'[\s_]+','', c.strip().lower()) for c in mraw.columns}
def find(name, *alts):
    for k,v in norm.items():
        if v==name or v in alts: return k
    return None
code_col    = find('code')
bs_col      = find('bodysite','bodysites','body_site')
type_col    = find('type','sampletype')
species_col = find('species')
dolphin_col = find('dolphinid','animalid','id')
meta = mraw.rename(columns={sid_col:'SampleID',
                            code_col:'code',
                            bs_col:'BodySite',
                            type_col:'Type',
                            species_col:'Species',
                            dolphin_col:'DolphinID'})
for c in ['SampleID','code','BodySite','Type','Species','DolphinID']:
    meta[c] = meta[c].astype(str).str.strip()

# Filters + pub labels
meta = meta[meta['BodySite'].str.lower()!='skin'].copy()
meta = meta[meta['DolphinID'].astype(str)!='57'].copy()
meta['BodySite_pub'] = normalize_bodysite_series(meta['BodySite'], meta['DolphinID'])
meta = meta[~meta['BodySite_pub'].isin({'Fetal wall','Maternal wall','Milk'})].copy()

def sp_norm(s):
    s = str(s).strip().upper()
    if s in {'TT','TURSIOPS','TURSIOPS TRUNCATUS'}: return 'Tt'
    if s in {'SB','STENO','STENO BREDANENSIS'}: return 'Sb'
    if s in {'TT/SB','MIXED'}: return 'Tt'
    return s.title() if s else s
meta['Species_norm']=meta['Species'].apply(sp_norm)

is_water    = meta['Type'].str.lower()=='water'
meta_nowater= meta[~is_water].copy()
meta_water  = meta[is_water].copy()

def code_norm_val(x):
    if x is None: return ''
    s=str(x).strip()
    if s=='' or s.lower() in {'na','n/a','n.a.','none','nan'}: return ''
    return re.sub(r'[^A-Za-z0-9]+','',s).upper()
meta_water = meta_water[meta_water['code'].apply(code_norm_val)!='']

# --- table ---
df = pd.read_csv(tsv_path, sep='\t', skiprows=1).rename(columns={'#OTU ID':'Phylum'})

def clean_phylum(s):
    s = str(s).strip()
    m = re.search(r'(?:^|[;|])p__\s*([^;|]+)', s)
    if m:
        ph = m.group(1).strip()
        if ph in {'', '_', '__', 'Unassigned'}:
            dom_m = re.search(r'(?:^|[;|])k__\s*([^;|]+)', s)
            dom = (dom_m.group(1).strip() if dom_m else '').lower()
            return 'Unclassified bacteria' if dom=='bacteria' else ('Unclassified archaea' if dom=='archaea' else 'Unclassified')
        return ph
    if ';' in s:
        parts = [p.strip().strip('_') for p in s.split(';')]
        dom = (parts[0] if parts else '').lower()
        phy = parts[1] if len(parts) > 1 else ''
        if phy == '' or phy == '_':
            if dom == 'bacteria': return 'Unclassified bacteria'
            if dom == 'archaea':  return 'Unclassified archaea'
            return 'Unclassified'
        return phy
    s2 = re.sub(r'^[A-Za-z]__', '', s).strip()
    return s2 if s2 else 'Unclassified'

df['Phylum'] = df['Phylum'].map(clean_phylum)
df = df.set_index('Phylum')

keep = set(meta['SampleID'])
df = df[[c for c in df.columns if c in keep]]

# --------------------------------------------------------------------
# Re-ordering + label/layout
# --------------------------------------------------------------------
codes_order = ['HI1','HI2','FL1','FL2','FL3','FL4']
BS_ORDER    = ['Vaginal','Penile','Mammary','Gastric','Rectal','Respiratory']

nonw = meta_nowater.copy()
is_fl4  = nonw['code'].str.upper().eq('FL4')
fl4_tt  = is_fl4 & nonw['Species_norm'].str.upper().eq('TT')
fl4_sb  = is_fl4 & nonw['Species_norm'].str.upper().eq('SB')

# Combined block
in_combined = (~is_fl4 & nonw['code'].str.upper().isin(codes_order[:-1])) | fl4_tt
combined_df = nonw[in_combined].copy()
combined_df.loc[fl4_tt[fl4_tt].index, 'code'] = 'FL4'

def bs_rank(bs):
    try: return BS_ORDER.index(bs)
    except ValueError: return 999

code_rank_map = {c:i for i,c in enumerate(codes_order)}
combined_df['_bsrank']   = combined_df['BodySite_pub'].map(bs_rank)
combined_df['_coderank'] = combined_df['code'].str.upper().map(code_rank_map)
combined_df = combined_df.sort_values(['_bsrank','_coderank','SampleID']).drop(columns=['_bsrank','_coderank'])

# Body-site runs for Combined
combined_runs=[]; cur=None; start=0
for i,bs in enumerate(combined_df['BodySite_pub'].tolist()):
    if cur is None: cur=bs; start=0
    elif bs!=cur: combined_runs.append((cur,start,i-1)); cur=bs; start=i
if len(combined_df)>0: combined_runs.append((cur,start,len(combined_df)-1))

# (BodySite, Code) contiguous segments inside Combined
segments = []
if not combined_df.empty:
    bs_prev = None; code_prev = None; a = 0
    for i,(bs,code) in enumerate(zip(combined_df['BodySite_pub'], combined_df['code'])):
        if i==0:
            bs_prev, code_prev = bs, code; continue
        if (bs != bs_prev) or (str(code).upper() != str(code_prev).upper()):
            segments.append((bs_prev, str(code_prev).upper(), a, i-1))
            a = i; bs_prev, code_prev = bs, code
    segments.append((bs_prev, str(code_prev).upper(), a, len(combined_df)-1))

# ---- CHANGED IN v3q: FL4_Sb block now respects BS_ORDER ----
fl4_sb_df = nonw[fl4_sb].copy()
fl4_sb_df['_bsrank'] = fl4_sb_df['BodySite_pub'].map(bs_rank)
fl4_sb_df = fl4_sb_df.sort_values(['_bsrank','SampleID']).drop(columns=['_bsrank'])

# Water blocks (unchanged)
def sort_by_pub(d): return d.sort_values(['BodySite_pub','SampleID'])
water_blocks=[]
for c in ['HI1','HI2','FL1','FL2','FL3','FL4']:
    part = meta_water[meta_water['code'].str.upper()==c].copy()
    if part.empty: continue
    part = sort_by_pub(part)
    water_blocks.append((f'Water_{c}', part))

# Assemble blocks
blocks = []
bs_runs = []

combined_sids = combined_df['SampleID'].tolist()
blocks.append(('Combined', combined_sids))
bs_runs.append(('Combined', combined_runs))

sb_sids = fl4_sb_df['SampleID'].tolist()
if sb_sids:
    blocks.append(('FL4_Sb', sb_sids))
    runs=[]; cur=None; s0=0
    for i,bs in enumerate(fl4_sb_df['BodySite_pub'].tolist()):
        if cur is None: cur=bs; s0=0
        elif bs!=cur: runs.append((cur,s0,i-1)); cur=bs; s0=i
    if len(fl4_sb_df)>0: runs.append((cur,s0,len(fl4_sb_df)-1))
    bs_runs.append(('FL4_Sb', runs))

for lbl,part in water_blocks:
    blocks.append((lbl, part['SampleID'].tolist()))
    runs=[]; cur=None; s0=0
    for i,bs in enumerate(part['BodySite_pub'].tolist()):
        if cur is None: cur=bs; s0=0
        elif bs!=cur: runs.append((cur,s0,i-1)); cur=bs; s0=i
    if len(part)>0: runs.append((cur,s0,len(part)-1))
    bs_runs.append((lbl, runs))

# --- x positions with gaps ---
gap    = 0.6
bs_gap = 0.6

x=[]; ranges=[]
pos=0.0
combined_start_index = None
for lbl,sids in blocks:
    n=len(sids)
    if n==0: continue
    s_index=len(x)
    if lbl=='Combined':
        combined_start_index = s_index
        prev_bs = None
        for sid, bs in zip(combined_sids, combined_df['BodySite_pub'].tolist()):
            if prev_bs is not None and bs != prev_bs:
                pos += bs_gap
            x.append(pos); pos += 1.0
            prev_bs = bs
        pos += gap
    else:
        for _ in range(n):
            x.append(pos); pos += 1.0
        pos += gap
    e_index=len(x)-1
    ranges.append((lbl, s_index, e_index))

x=np.array(x)

ordered=[sid for _,lst in blocks for sid in lst]
df = df[ordered]

TOP=20
totals=df.sum(axis=1).sort_values(ascending=False)
plot_df=df.loc[totals.head(TOP).index].copy()
other=(df.sum(axis=0)-plot_df.sum(axis=0)).clip(lower=0)
plot_df.loc['Other phyla']=other
plot_df=plot_df.loc[list(totals.head(TOP).index)+['Other phyla']]

pal=list(plt.get_cmap('tab20').colors)
while len(pal)<plot_df.shape[0]: pal+=list(plt.get_cmap('tab20b').colors)
colors={p:pal[i%len(pal)] for i,p in enumerate(plot_df.index)}

fig,ax=plt.subplots(figsize=(22,8))
fig.subplots_adjust(top=0.66)

bottom=np.zeros(len(x))
for phy in plot_df.index[::-1]:
    vals=plot_df.loc[phy,ordered].values
    ax.bar(x,vals,width=0.9,bottom=bottom,color=colors[phy],linewidth=0)
    bottom+=vals

trans = blended_transform_factory(ax.transData, ax.transAxes)
yl_top     = 1.24
yl_bottom  = 1.04
yl_group   = 1.34
FS_BODY    = 12
FS_CODE    = 11
FS_HEADER  = 12

if combined_start_index is not None:
    for bs, code, a_loc, b_loc in segments:
        s = combined_start_index + a_loc
        e = combined_start_index + b_loc
        ax.text((x[s]+x[e])/2.0, yl_bottom, code, transform=trans,
                rotation=90, ha='center', va='bottom', fontsize=FS_CODE,
                fontweight='normal', clip_on=False, zorder=6)

for lbl,s,e in ranges:
    if lbl.startswith('Water_'):
        code = lbl.split('_',1)[1]
        ax.text((x[s]+x[e])/2.0, yl_bottom, code, transform=trans,
                rotation=90, ha='center', va='bottom', fontsize=FS_CODE,
                fontweight='normal', clip_on=False, zorder=5)

offset={}; cur=0
for lbl,sids in blocks:
    offset[lbl]=cur; cur+=len(sids)

for lbl,runs in bs_runs:
    base=offset.get(lbl,0)
    for (bs,a,b) in runs:
        if a>b: continue
        s = base+a; e = base+b
        cx = (x[s]+x[e])/2.0
        if lbl=='Combined':
            ax.text(cx, yl_top, bs, transform=trans,
                    ha='center', va='bottom', fontsize=FS_BODY, fontweight='bold',
                    rotation=0, clip_on=False, zorder=6)
        elif lbl=='FL4_Sb':
            ax.text(cx, yl_bottom, bs, transform=trans,
                    ha='center', va='bottom', fontsize=FS_BODY-1, fontweight='normal',
                    rotation=90, clip_on=False, zorder=6)
        else:
            ax.text(cx, yl_top, bs, transform=trans,
                    ha='center', va='bottom', fontsize=FS_BODY, fontweight='normal',
                    rotation=90, clip_on=False, zorder=6)

water_spans = [(s,e) for (lbl,s,e) in ranges if lbl.startswith('Water_')]
if water_spans:
    s_all = min(s for s,_ in water_spans)
    e_all = max(e for _,e in water_spans)
    ax.text((x[s_all]+x[e_all])/2.0, yl_top, 'Water', transform=trans,
            ha='center', va='bottom', fontsize=FS_HEADER, fontweight='bold',
            clip_on=False, zorder=7)

for lbl,s,e in ranges:
    if lbl=='FL4_Sb':
        ax.text((x[s]+x[e])/2.0, yl_group, 'Sb', transform=trans,
                ha='center', va='bottom', fontsize=FS_HEADER, fontweight='bold',
                clip_on=False, zorder=7)
        break

ax.set_ylim(0,1.0)
ax.set_ylabel('Relative abundance')
ax.set_xlim(x.min()-0.6, x.max()+0.6)
ax.set_xticks([])
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(False)

handles=[plt.Rectangle((0,0),1,1,color=colors[p]) for p in plot_df.index]
ax.legend(handles, list(plot_df.index), bbox_to_anchor=(1.01,1), loc='upper left',
          frameon=False, fontsize=9)

# ---- Save figures under v3q name ----
for ext in ('png','pdf','svg'):
    plt.savefig(os.path.join(OUT, f'phylum_reorder_by_code_no57_v3q.{ext}'),
                dpi=300 if ext=='png' else None, bbox_inches='tight')

# ---- Export the sample order (v3q) ----
pd.DataFrame({'SampleID': [sid for _,lst in blocks for sid in lst]}).to_csv(
    os.path.join(OUT,'sample_order_reorder_by_code_no57_v3q.tsv'), sep='\t', index=False)

from matplotlib.text import Text
fig_for_report = plt.gcf()
def export_font_report(fig, report_path, extra=None):
    lines = []
    lines.append(f"Matplotlib version: {matplotlib.__version__}")
    lines.append(f"rcParams[font.family]: {mpl.rcParams.get('font.family')}")
    lines.append(f"rcParams[font.sans-serif]: {mpl.rcParams.get('font.sans-serif')}")
    lines.append(f"User-registered Arial font files: {len(FOUND)}")
    for p in FOUND:
        lines.append(f"  - {p}")
    resolved = findfont(FontProperties(family=['Arial','sans-serif']))
    lines.append(f"Resolved font path for 'Arial' request: {resolved}\n")

    names = []
    for i, t in enumerate(fig.findobj(match=Text), 1):
        fp = t.get_fontproperties()
        name = fp.get_name()
        size = fp.get_size_in_points()
        sample = (t.get_text() or "").replace("\n"," ")[:60]
        names.append(name)
        lines.append(f"[Text#{i:03d}] name='{name}', size={size:.1f}pt, sample='{sample}'")

    lines.append("\n=== Unique font names & counts ===")
    for n, c in Counter(names).most_common():
        lines.append(f"{n}: {c}")

    if extra:
        lines.append("\n=== Extra info ===")
        for k, v in extra.items():
            lines.append(f"{k}: {v}")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

report_path = os.path.join(OUT, "font_report_phylum_v3q.txt")
export_font_report(fig_for_report, report_path, {"output_png": os.path.join(OUT, "phylum_reorder_by_code_no57_v3q.png")})
print("Wrote font report:", report_path)

PY

echo "Done."
