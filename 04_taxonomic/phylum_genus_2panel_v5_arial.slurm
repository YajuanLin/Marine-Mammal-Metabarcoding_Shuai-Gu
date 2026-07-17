#!/bin/bash
#SBATCH --job-name=phylum_genus_2panel_v5_arial
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/barplot/pick/phylum_genus_2panel_v5_arial.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/barplot/pick/phylum_genus_2panel_v5_arial.err
#SBATCH --time=00:35:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
OUT=$BASE/barplot/pick
INP=$BASE/barplot/all
META=$BASE/merged_dolphin_4.txt

PHYLUM_RELQZA_IN=$INP/phylum-relfreq.qza
PHYLUM_TABLE_QZA=$INP/phylum-table.qza
PHYLUM_TSV=$OUT/phylum_relab_raw.tsv

GENUS_RELQZA_IN=$BASE/genus-relfreq.qza
GENUS_TABLE_QZA=$BASE/genus-table.qza
GENUS_TSV=$OUT/genus_relab.tsv

mkdir -p "$OUT"

# ---- phylum tsv (same logic as v16) ----
if [ ! -f "$PHYLUM_TSV" ]; then
  if [ -f "$PHYLUM_RELQZA_IN" ]; then
    rm -rf "$OUT/export_phylum"
    qiime tools export --input-path "$PHYLUM_RELQZA_IN" --output-path "$OUT/export_phylum"
    biom convert -i "$OUT/export_phylum/feature-table.biom" -o "$PHYLUM_TSV" --to-tsv
  elif [ -f "$PHYLUM_TABLE_QZA" ]; then
    RELQZA_TMP="$OUT/phylum-relfreq.tmp.qza"
    qiime feature-table relative-frequency --i-table "$PHYLUM_TABLE_QZA" --o-relative-frequency-table "$RELQZA_TMP"
    rm -rf "$OUT/export_phylum"
    qiime tools export --input-path "$RELQZA_TMP" --output-path "$OUT/export_phylum"
    biom convert -i "$OUT/export_phylum/feature-table.biom" -o "$PHYLUM_TSV" --to-tsv
  else
    echo "Missing phylum inputs: $PHYLUM_RELQZA_IN or $PHYLUM_TABLE_QZA" >&2
    exit 1
  fi
fi

# ---- genus tsv (same logic as v15legend6) ----
if [ ! -f "$GENUS_TSV" ]; then
  if [ ! -f "$GENUS_RELQZA_IN" ]; then
    qiime feature-table relative-frequency \
      --i-table "$GENUS_TABLE_QZA" \
      --o-relative-frequency-table "$GENUS_RELQZA_IN"
  fi
  rm -rf "$OUT/export_genus"
  qiime tools export --input-path "$GENUS_RELQZA_IN" --output-path "$OUT/export_genus"
  biom convert -i "$OUT/export_genus/feature-table.biom" -o "$GENUS_TSV" --to-tsv
fi

python - <<'PY'
import os, re, numpy as np, pandas as pd, matplotlib, shutil
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory as blend
from matplotlib import patheffects as pe
from matplotlib import font_manager as fm
from matplotlib import rcParams
from collections import defaultdict, Counter

# ---------- Arial font setup (minimal + persistent) ----------
# 1) Clear cache (safe if not present)
try:
    shutil.rmtree(os.path.expanduser('~/.cache/matplotlib'), ignore_errors=True)
except Exception:
    pass

# 2) Register local Arial TTFs
BASE = os.path.expanduser('~/sg/formal/merge4tubes')
OUT  = os.path.join(BASE,'barplot','pick')
arial_dir = os.path.join(BASE, 'fonts', 'arial')
arial_files = [
    'Arial.ttf',
    'Arial Bold.ttf',
    'Arial Italic.ttf',
    'Arial Bold Italic.ttf'
]
registered = []
for fn in arial_files:
    fp = os.path.join(arial_dir, fn)
    if os.path.exists(fp):
        try:
            fm.fontManager.addfont(fp)
            registered.append(fp)
        except Exception:
            pass

# 3) Global rcParams
rcParams['font.family'] = ['sans-serif']
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans', 'Noto Sans', 'Nimbus Sans', 'sans-serif']
# Mathtext: there is no true Arial math; dejavusans keeps appearance consistent and avoids Type 3
rcParams['mathtext.fontset'] = 'dejavusans'
rcParams['pdf.fonttype'] = 42       # embed Truetype
rcParams['svg.fonttype'] = 'none'   # keep text as text

import matplotlib
import matplotlib.font_manager as _fm

# ---------- paths & IO ----------
meta_path    = os.path.join(BASE,'merged_dolphin_4.txt')
phylum_tsv   = os.path.join(OUT,'phylum_relab_raw.tsv')
genus_tsv    = os.path.join(OUT,'genus_relab.tsv')

# ===== metadata (shared) =====
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
                            code_col or 'code':'code',
                            bs_col or 'BodySite':'BodySite',
                            type_col or 'Type':'Type',
                            species_col or 'Species':'Species',
                            dolphin_col or 'DolphinID':'DolphinID'}).copy()
for c in ['SampleID','code','BodySite','Type','Species','DolphinID']:
    meta[c] = meta[c].astype(str).str.strip()

def ci(s): return s.casefold() if isinstance(s,str) else s
def pick_by_type_blank(mt):
    return mt[mt['Type'].map(lambda x: ci(x)=='blank')].sort_values('SampleID')
def pick_by_bs(mt, name):
    return mt[mt['BodySite'].map(lambda x: ci(x)==ci(name))].sort_values('SampleID')
def norm_bs_for_order(x):
    xci = str(x).casefold()
    if xci=='gastric (mouth)': return 'Gastric'
    if xci=='genital p':       return 'Genital_P'
    if xci=='genital v':       return 'Genital_V'
    if xci=='blowhole':        return 'Respiratory'
    return str(x)
def pick_by_dolphin(mt, did, order_sites):
    d = mt[mt['DolphinID'].astype(str)==str(did)].copy()
    d['__bs__'] = d['BodySite'].map(norm_bs_for_order)
    out=[]
    for s in order_sites:
        alias = 'Respiratory' if ci(s) in ('blowhole','respiratory') else s
        part = d[d['__bs__']==alias].sort_values('SampleID')
        out.append(part)
    res = pd.concat(out, axis=0) if out else d.iloc[0:0]
    return res.drop(columns='__bs__', errors='ignore')
def pick_by_sampleid(mt, sid):
    return mt[mt['SampleID']==sid].copy()

# shared blocks/order
blocks=[]
def add_block(label, dfm):
    if not dfm.empty: blocks.append((label, dfm['SampleID'].tolist()))
add_block('Blank', pick_by_type_blank(meta))
add_block('Water', pick_by_sampleid(meta,'P8'))
add_block('External Placenta', pick_by_bs(meta,'External Placenta'))
add_block('Internal Placenta', pick_by_bs(meta,'Internal Placenta'))
add_block('Mother',  pick_by_dolphin(meta,11,['Genital_V','Mammary','Gastric','Rectal','Respiratory']))
add_block('Milk',    pick_by_bs(meta,'Milk'))
add_block('Neonate', pick_by_dolphin(meta,57,['Genital_V','Gastric','Respiratory']))
ordered = [sid for _, sids in blocks for sid in sids]

# ===== phylum data =====
df_phy = pd.read_csv(phylum_tsv, sep='\t', skiprows=1).rename(columns={'#OTU ID':'Phylum'})
def clean_phylum(x: str) -> str:
    s = str(x).strip()
    m = re.search(r'p__([^;|]+)', s)
    if m:
        name = m.group(1).strip()
        return 'Unassigned bacteria' if name in {'', '_', '__', 'unassigned'} else name
    parts = [t.strip() for t in re.split(r'[;|]', s) if t.strip()!='']
    if parts:
        last = parts[-1]
        if last.lower() in {'_', '__', 'unassigned', 'unclassified', 'uncultured'}:
            return 'Unassigned bacteria' if any('bacteria' in p.lower() for p in parts) else 'Unassigned'
        return last
    return s
df_phy['Phylum']=df_phy['Phylum'].map(clean_phylum)
df_phy = df_phy.set_index('Phylum')
df_phy = df_phy[[c for c in df_phy.columns if c in ordered]]

# ===== genus data =====
raw_g = pd.read_csv(genus_tsv, sep='\t', skiprows=1).rename(columns={'#OTU ID':'Taxon'}).fillna('')
def clean_genus(x):
    x=str(x)
    m=re.search(r'g__([^;|]+)', x)
    if m:
        g=m.group(1).strip()
        return g if g else 'Unassigned'
    x=re.sub(r'^[a-z]__','', x).strip()
    return x if x else 'Unassigned'
raw_g['Genus'] = raw_g['Taxon'].map(clean_genus)
sample_cols = [c for c in raw_g.columns if c in ordered]
df_gen = raw_g[['Genus'] + sample_cols].copy().set_index('Genus')

# ===== tops + Other =====
TOP_PHY, TOP_GEN = 20, 30
tot_phy = df_phy.sum(axis=1).sort_values(ascending=False)
plot_phy = df_phy.loc[tot_phy.head(TOP_PHY).index].copy()
other_phy = (df_phy.sum(axis=0) - plot_phy.sum(axis=0)).clip(lower=0)
plot_phy.loc['Other phyla'] = other_phy
plot_phy = plot_phy.loc[list(tot_phy.head(TOP_PHY).index)+['Other phyla']]

tot_gen = df_gen.sum(axis=1).sort_values(ascending=False)
plot_gen = df_gen.loc[tot_gen.head(TOP_GEN).index].copy()
other_gen = (df_gen.sum(axis=0) - plot_gen.sum(axis=0)).clip(lower=0)
plot_gen.loc['Other genera'] = other_gen
plot_gen = plot_gen.loc[list(tot_gen.head(TOP_GEN).index)+['Other genera']]

# ===== shared x positions =====
gap=0.6
x=[]; ranges=[]; start=0.0
for lbl,sids in blocks:
    n=len(sids)
    if n>0:
        xs=[start+j for j in range(n)]
        ranges.append((lbl, len(x), len(x)+n-1))
        x += xs
        start += n + gap
x = np.array(x, dtype=float)
def center(s,e): return (x[s]+x[e])/2.0

# ===== colors =====
def make_colors(index):
    pal=list(plt.get_cmap('tab20').colors)
    while len(pal) < len(index): pal += list(plt.get_cmap('tab20b').colors)
    return {k: pal[i % len(pal)] for i,k in enumerate(index)}
colors_phy = make_colors(plot_phy.index)
colors_gen = make_colors(plot_gen.index)

# ===== layout (two rows + legend columns) =====
fig = plt.figure(figsize=(22, 15))
gs = plt.GridSpec(2, 2, figure=fig, width_ratios=[1.0, 0.42], height_ratios=[1.0, 1.0], wspace=0.05, hspace=0.25)
ax_phy   = fig.add_subplot(gs[0,0])
ax_phy_L = fig.add_subplot(gs[0,1]); ax_phy_L.set_axis_off()
ax_gen   = fig.add_subplot(gs[1,0])
ax_gen_L = fig.add_subplot(gs[1,1]); ax_gen_L.set_axis_off()

# ===== plot helper =====
def plot_stacked(ax, plot_df, colors):
    bottom = np.zeros(len(x), dtype=float)
    for key in list(plot_df.index)[::-1]:
        vals = plot_df.loc[key, ordered].values.astype(float)
        ax.bar(x, vals, width=0.9, bottom=bottom, color=colors[key], linewidth=0)
        bottom += vals
    ax.set_ylim(0,1.0); ax.grid(False)
    ax.set_xlim(x.min()-0.6, x.max()+0.6)
    ax.set_xticks([])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_ylabel('Relative abundance')

# headers/subtitles
data_top_phy = blend(ax_phy.transData, ax_phy.transAxes)
data_top_gen = blend(ax_gen.transData, ax_gen.transAxes)
yl_hi, yl_lo = 1.12, 1.10
min_dx = 1.4
def wrap_label(lbl: str) -> str:
    if lbl == 'External Placenta': return 'Maternal\nwall\n(Placenta)'
    if lbl == 'Internal Placenta': return 'Fetal\nwall\n(Placenta)'
    return lbl
def draw_headers(ax, data_top):
    centers = [(center(s,e), lbl) for (lbl,s,e) in ranges]
    centers.sort(key=lambda t: t[0])
    row_toggle = 0
    for i,(cx,lbl) in enumerate(centers):
        if i>0 and abs(cx - centers[i-1][0]) < min_dx:
            row_toggle = 1 - row_toggle
        else:
            row_toggle = 0
        y = yl_hi if row_toggle else yl_lo
        ax.text(cx, y, wrap_label(lbl), ha='center', va='bottom',
                fontsize=12, fontweight='bold', rotation=0, transform=data_top)
def bs_display(name: str, *, keep_blowhole: bool = False) -> str:
    if not isinstance(name, str): return name
    n = name.casefold().strip()
    if n in ('genital_v','genital v'): return 'Vaginal'
    if n in ('genital_p','genital p'): return 'Penile'
    if n == 'blowhole': return 'Blowhole' if keep_blowhole else 'Respiratory'
    if n == 'respiratory': return 'Respiratory'
    if n == 'gastric (mouth)': return 'Oral'
    if n == 'gastric': return 'Gastric'
    if n == 'external placenta': return 'Maternal wall'
    if n == 'internal placenta': return 'Fetal wall'
    return name
def draw_subtitles(ax, data_top):
    subtitle_y = 1.045
    meta_ix = meta.set_index('SampleID')
    for lbl, sids in blocks:
        if lbl in {'Mother','Neonate'} and len(sids) > 0:
            sub = meta_ix.loc[sids, 'BodySite'].astype(str).tolist()
            dids = meta_ix.loc[sids, 'DolphinID'].astype(str).tolist()
            keep_blowhole = any(d == '57' for d in dids)
            disp = [bs_display(b, keep_blowhole=keep_blowhole) for b in sub]
            starts=[0]
            for i in range(1,len(disp)):
                if disp[i]!=disp[i-1]: starts.append(i)
            starts.append(len(disp))
            baseS, baseE = next((S,E) for (L,S,E) in ranges if L==lbl)
            for k in range(1,len(starts)):
                a=starts[k-1]; b=starts[k]-1
                xa = x[baseS + a]; xb = x[baseS + b]
                ax.text((xa+xb)/2.0, subtitle_y, disp[a], ha='center', va='bottom',
                        fontsize=9, rotation=0, transform=data_top)

# draw both panels
plot_stacked(ax_phy, plot_phy, colors_phy)
draw_headers(ax_phy, data_top_phy)
draw_subtitles(ax_phy, data_top_phy)

plot_stacked(ax_gen, plot_gen, colors_gen)
draw_headers(ax_gen, data_top_gen)
draw_subtitles(ax_gen, data_top_gen)

# ===== Panel labels outside plotting area (above each y-axis) =====
fig = plt.gcf()
fig.canvas.draw()
offset_y = 0.008
offset_x = 0.002
for ax, letter in [(ax_phy, 'A'), (ax_gen, 'B')]:
    bb = ax.get_position(fig)
    tx = fig.text(bb.x0 + offset_x, bb.y1 + offset_y, letter,
                  ha='left', va='bottom', fontsize=16, fontweight='bold')
    tx.set_path_effects([pe.withStroke(linewidth=3, foreground='white')])

# legends
handles_phy = [plt.Rectangle((0,0),1,1, facecolor=colors_phy[k]) for k in plot_phy.index]
ax_phy_L.legend(handles_phy, list(plot_phy.index), loc='upper left', frameon=False, fontsize=9)

rows = raw_g[['Genus','Taxon'] + sample_cols].copy()
rows['_sum'] = rows[sample_cols].sum(axis=1).astype(float)
rows = rows.sort_values('_sum', ascending=False)
tax_by_genus = defaultdict(list)
for _, r in rows.iterrows():
    tax_by_genus[r['Genus']].append(r['Taxon'])
def _parse_tax(line: str):
    s=(line or '').strip()
    ranks={"k":"","p":"","c":"","o":"","f":"","g":""}
    found=False
    for r in ranks:
        m=re.search(rf"(?:^|[;|]){r}__\s*([^;|]+)",s)
        if m: ranks[r]=m.group(1).strip().strip('_'); found=True
    if not found:
        parts=[p.strip().strip('_') for p in s.split(';') if p.strip()]
        if parts:
            ranks["k"]=parts[0]
            for i,r in enumerate(["p","c","o","f","g"],start=1):
                if len(parts)>i: ranks[r]=parts[i]
    return ranks
def _nearest_non_is_parent(r):
    for key in ["f","o","c","p","k"]:
        v=(r.get(key,"") or "").strip()
        if v and v.lower()!="incertae sedis" and v not in {"_","__","Unassigned"}:
            return v
    return "Bacteria"
def legend_label_from_taxon(tax_line: str):
    r=_parse_tax(tax_line)
    g=(r.get("g","") or "").strip()
    if g and g.lower() not in {"incertae sedis","unassigned","__","_"}:
        return g
    parent=_nearest_non_is_parent(r)
    if g.lower()=="incertae sedis":
        return f"Incertae sedis (within {parent})"
    return f"Unassigned {parent}"
consume=defaultdict(int); labels=[]
for key in list(plot_gen.index):
    if key=="Other genera":
        labels.append("Other genera"); continue
    lst=tax_by_genus.get(key, [])
    k=consume[key]; consume[key]+=1
    tax=lst[k] if k < len(lst) else ""
    labels.append(legend_label_from_taxon(tax))
from collections import Counter
cnt=Counter(l for l in labels if l.startswith("Incertae sedis"))
seen=defaultdict(int)
labels=[ (f"{l} [{(seen.__setitem__(l, seen[l]+1) or seen[l])}]" if l.startswith("Incertae sedis") and cnt[l]>1 else l)
         for l in labels ]
handles_gen=[plt.Rectangle((0,0),1,1, facecolor=colors_gen[k]) for k in plot_gen.index]
ax_gen_L.legend(handles_gen, labels, loc='upper left', frameon=False, fontsize=9)

# save
for ext in ('png','pdf','svg'):
    fig.savefig(os.path.join(OUT, f'phylum_genus_2panel_v5_arial.{ext}'),
                dpi=300 if ext=='png' else None, bbox_inches='tight')
print("Wrote:", os.path.join(OUT, 'phylum_genus_2panel_v5_arial.*'))

# ---------- font report ----------
try:
    report = os.path.join(OUT, 'phylum_genus_2panel_v5_arial_font_report.txt')
    with open(report, 'w') as f:
        f.write(f"Matplotlib version: {matplotlib.__version__}\n")
        f.write(f"rcParams[font.family]: {rcParams.get('font.family')}\n")
        f.write(f"rcParams[font.sans-serif]: {rcParams.get('font.sans-serif')}\n")
        f.write(f"mathtext.fontset = {rcParams.get('mathtext.fontset')}\n")
        f.write(f"pdf.fonttype = {rcParams.get('pdf.fonttype')}\n")
        f.write(f"svg.fonttype = {rcParams.get('svg.fonttype')}\n")
        f.write(f"User-registered Arial font files: {len(registered)}\n")
        for p in registered:
            f.write(f"  - {p}\n")
        resolved = _fm.findfont('Arial')
        f.write(f"Resolved font path for 'Arial' request: {resolved}\n")
except Exception:
    pass
PY

echo "Done."
