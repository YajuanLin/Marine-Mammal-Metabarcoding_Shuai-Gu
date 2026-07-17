#!/bin/bash
#SBATCH --job-name=genus-heatmap-onepanel-v3p
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/heatmap/genus/genus_heatmap_onepanel_v3p.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/heatmap/genus/genus_heatmap_onepanel_v3p.err
#SBATCH --time=01:20:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu
set -euo pipefail
module load qiime2/2024.2

BASE=~/sg/formal/merge4tubes
OUTD=$BASE/heatmap/genus
mkdir -p "$OUTD"

# clear per-user Matplotlib cache to avoid stale fallbacks
rm -rf ~/.cache/matplotlib || true

python - <<'PY'
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib as mpl
from matplotlib import font_manager as fm
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.transforms import blended_transform_factory as blend
from collections import Counter

# ================= Fonts: force ARIAL everywhere =================
# Register user-provided Arial fonts (you uploaded these already)
FONT_DIR = os.path.expanduser("~/sg/formal/merge4tubes/fonts/arial")
FOUND = []
if os.path.isdir(FONT_DIR):
    for fn in os.listdir(FONT_DIR):
        p = os.path.join(FONT_DIR, fn)
        if p.lower().endswith(".ttf") and "arial" in fn.lower():
            try:
                fm.fontManager.addfont(p)
                FOUND.append(p)
            except Exception:
                pass
# Rebuild cache so new fonts are visible
try:
    fm._rebuild()
except Exception:
    pass

# Set global rcParams so ALL text uses Arial (including mathtext)
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial","Helvetica","Liberation Sans","DejaVu Sans","Noto Sans","Nimbus Sans","sans-serif"],
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "pdf.fonttype": 42,   # embed TrueType in PDF
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})

# ================= Config =================
BASE = os.path.expanduser('~/sg/formal/merge4tubes')
META = os.path.join(BASE,'merged_dolphin_4.txt')
OUTD = os.path.join(BASE,'heatmap/genus')
os.makedirs(OUTD, exist_ok=True)

bodysite_order = ['Vaginal','Penile','Mammary','Gastric','Rectal','Respiratory']
slug = {'Vaginal':'vaginal','Penile':'penile','Mammary':'mammary',
        'Gastric':'gastric','Rectal':'rectal','Respiratory':'respiratory'}
species_order  = ['Tt','Sb']
facility_order = ['HI1','HI2','FL1','FL2','FL3','FL4']
fac_colors = {'HI1':'#4477AA','HI2':'#EE6677','FL1':'#228833','FL2':'#AA3377','FL3':'#66CCEE','FL4':'#CCBB44','Unknown':'#999999'}

SPACER_COLS_SPECIES = 1
SPACER_COLS_BODYSITE = 4

dpi, fig_w = 300, 18
row_height = 0.24
fs_genus = 11
fs_site  = 15
fw_site  = 'bold'

# ---- LOG10 colormap ----
PSEUDO_LOG = 1e-4
cmap_wr = LinearSegmentedColormap.from_list('whitered', ['#FFFFFF', '#B2182B'])

# ================= Utilities =================
def pick(df,t,*alts):
    n={c:re.sub(r'[\s_]+','',c.strip().lower()) for c in df.columns}
    for k,v in n.items():
        if v==t or v in alts: return k
    return None

def norm_sp(x):
    s=str(x).strip().lower()
    if s.startswith('t'): return 'Tt'
    if s.startswith('s'): return 'Sb'
    return s

def bs_norm(raw):
    s=re.sub(r'[\s_]+','',str(raw).strip().lower())
    if s in {'genitalv','genital-v','vaginal'}: return 'Vaginal'
    if s in {'genitalp','genital-p','penile'}:  return 'Penile'
    if 'mammary' in s:                           return 'Mammary'
    if 'gastric'  in s:                           return 'Gastric'
    if 'rectal'   in s:                           return 'Rectal'
    if ('blowhole' in s) or ('respir' in s):      return 'Respiratory'
    return ''

def normalize_semicolon(t):
    s = str(t).strip()
    parts = [p.strip() for p in s.split(';')]
    if len(parts) < 6: parts += [''] * (6 - len(parts))
    parts = parts[:6]
    parts = [p if p not in {'','_','__'} else '__' for p in parts]
    return ';'.join(parts)

def split_semicolon(s):
    parts = str(s).split(';')
    parts = (parts + ['__'] * 6)[:6]
    return parts  # K,P,C,O,F,G

def _first_parent_nonincertae(*levels):
    for x in levels:
        if x and x != '__' and 'incertae' not in x.lower():
            return x
    for x in levels:
        if x and x != '__':
            return x
    return None

def pretty_from_semicolon(s):
    K,P,C,O,F,G = split_semicolon(s)
    g = (G or '').strip()
    if g and g != '__' and 'incertae' not in g.lower():
        return g
    parent = _first_parent_nonincertae(F, O, C, P, K)
    if g and 'incertae' in g.lower():
        return f'Incertae sedis (within {parent})' if parent else 'Incertae sedis'
    return f'Unassigned {parent}' if parent else 'Unassigned'

GENUS_EXCEPTIONS_PREFIXES = ('Incertae sedis','Unassigned','UCG','NS','Clade','Cluster','Group')

def _is_family_or_higher(word: str) -> bool:
    return bool(re.search(r'(aceae|ales|bacteria|archaea|mycetes|phyta|opsida|phyceae)$', word))

def split_genus_parenthetical(label: str):
    s = label.strip()
    if s.startswith("Incertae sedis"):
        m = re.match(r"^(Incertae sedis)(.*)$", s)
        pre = m.group(1); suf = m.group(2) or ""
        return (pre, "", suf, 'incertae')
    if any(s.startswith(p) for p in GENUS_EXCEPTIONS_PREFIXES if p != 'Incertae sedis'):
        return ("", "", s, 'plain')
    m = re.match(r"^Candidatus(\s+)([A-Z][a-z]+)(.*)$", s)
    if m:
        space = m.group(1); genus = m.group(2); rest = m.group(3)
        return ("Candidatus"+space, genus, rest, 'cand')
    m = re.match(r"^([A-Z][a-z]+(?:-[A-Z][a-z]+)+)(\s*(\(.*\))?)$", s)
    if m:
        return ("", m.group(1), m.group(2) or '', 'genus')
    m = re.match(r"^([A-Z][a-z]+)(\s*(\(.*\))?)$", s)
    if m:
        genus = m.group(1); rest = m.group(2) or ''
        if not _is_family_or_higher(genus):
            return ("", genus, rest, 'genus')
    return ("", "", s, 'plain')

# ================= Load metadata =================
m = pd.read_csv(META, sep='\t', dtype=str).fillna('')
sid=m.columns[0]; code=pick(m,'code'); bs=pick(m,'bodysite','bodysites','body_site')
sp=pick(m,'species','sp'); did=pick(m,'dolphinid','animalid','id'); typ=pick(m,'type','sampletype')
m = m.rename(columns={sid:'SampleID',code:'code',bs:'BodySite',sp:'Species',did:'DolphinID',typ:'Type'})
for c in ['SampleID','code','BodySite','Species','DolphinID','Type']: m[c]=m[c].astype(str).str.strip()
m['SP2']=m['Species'].map(norm_sp)
m['BS2']=m['BodySite'].map(bs_norm)
m = m[(m['DolphinID']!='57') & (m['Type'].str.lower()!='water') &
      (m['SP2'].isin(['Tt','Sb'])) & (m['BS2']!='')].copy()

# ================= Genus union from barplots =================
wanted, seen, missing = [], set(), []
for bs in bodysite_order:
    f = os.path.join(BASE,'barplot',slug[bs], f'genus_drawn_{slug[bs]}_v13.tsv')
    if not os.path.isfile(f): missing.append(f); continue
    raw = pd.read_csv(f, sep='\t', dtype=str)['Genus'].tolist()
    for x in raw:
        if x == 'Other genera': continue
        g = normalize_semicolon(x)
        if g not in seen:
            seen.add(g); wanted.append(g)
if not wanted:
    raise SystemExit("No genera from genus_drawn_*_v13.tsv after normalization.")
if missing:
    print("WARNING missing genus_drawn files:", *missing, sep='\n- ')

# ================= Per-site matrices =================
site_mats = {}
for bs in bodysite_order:
    tsv = os.path.join(BASE,'barplot',slug[bs],'genus_relab.tsv')
    if not os.path.isfile(tsv):
        raise SystemExit(f"Missing TSV: {tsv}")
    df = pd.read_csv(tsv, sep='\t', skiprows=1).rename(columns={'#OTU ID':'Tax'})
    df['Genus'] = df['Tax'].map(normalize_semicolon)
    df = df.drop(columns=['Tax']).set_index('Genus').apply(pd.to_numeric, errors='coerce').fillna(0.0)
    df = df.groupby(level=0).sum()
    keep = set(m.loc[m['BS2']==bs,'SampleID'])
    df = df[[c for c in df.columns if c in keep]]
    site_mats[bs] = df

# ================= Column order =================
ordered_cols, block_bounds = [], []
for bs in bodysite_order:
    sub = m[m['BS2']==bs]
    bs_cols=[]
    for spc in species_order:
        sub_sp=sub[sub['SP2']==spc]
        for fac in facility_order:
            part=sub_sp[sub_sp['code']==fac].sort_values('SampleID')
            if not part.empty: bs_cols += part['SampleID'].tolist()
        if spc != species_order[-1] and len(sub_sp)>0:
            bs_cols += [None]*SPACER_COLS_SPECIES
    start=len(ordered_cols); ordered_cols+=bs_cols; end=len(ordered_cols)-1
    block_bounds.append((start,end,bs))
    if bs != bodysite_order[-1]: ordered_cols += [None]*SPACER_COLS_BODYSITE
real_samples = [c for c in ordered_cols if c is not None]

# ================= Global matrix =================
X = pd.DataFrame(0.0, index=wanted, columns=real_samples)
for bs in bodysite_order:
    df = site_mats[bs].reindex(index=wanted, fill_value=0.0)
    cols_here = [c for c in real_samples if c in df.columns]
    if cols_here:
        X.loc[:, cols_here] = df[cols_here].values

# ================= Labels & highlight set =================
def pretty_from_semicolon_cached(s):
    return pretty_from_semicolon(s)

pretty_labels = [pretty_from_semicolon_cached(s) for s in X.index]
HIGHLIGHT_GENUS_INPUT = [
    "Aggregatibacter","Haemophilus","Staphylococcus","Fusobacterium",
    "Leptotrichia","Mycoplasmoides","Arcobacter","Helicobacter",
    "Salinivibrio","Incertae sedis (within Spirochaetaceae)","Ligilactobacillus",
]
BOLD_SET   = {s.strip().lower() for s in HIGHLIGHT_GENUS_INPUT}
GRAY_COLOR = "#7a7a7a"

# ================= Transform & draw =================
PSEUDO_LOG = 1e-4
A = X.to_numpy()
T = np.log10(A + PSEUDO_LOG)
vmin, vmax = np.log10(PSEUDO_LOG), 0.0
cmap_wr = LinearSegmentedColormap.from_list('whitered', ['#FFFFFF', '#B2182B'])
cmap = cmap_wr.copy(); cmap.set_bad('white')

# Matrix with NaN spacers
mat_cols=[]
for sid in ordered_cols:
    mat_cols.append(np.full((T.shape[0],), np.nan) if sid is None else T[:, real_samples.index(sid)])
H = np.column_stack(mat_cols)
Hm = np.ma.masked_invalid(H)

# Facility strip
def hex2rgb01(h): return tuple(int(h[i:i+2],16)/255 for i in (1,3,5))
sid_to_fac = dict(zip(m['SampleID'], m['code']))
fac_row = [(1,1,1) if sid is None else hex2rgb01(fac_colors.get(sid_to_fac.get(sid,'Unknown'),'#999999'))
           for sid in ordered_cols]
fac_strip = np.array([fac_row])

# Layout
dpi, fig_w = 300, 18
row_height = 0.24
fs_genus = 11
fs_site  = 15
fw_site  = 'bold'

fig_h = max(6, row_height * X.shape[0])
fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)

ax  = fig.add_axes([0.19, 0.12, 0.78, 0.78])
cax = fig.add_axes([0.105, 0.22, 0.018, 0.26])
leg_ax = fig.add_axes([0.105, 0.56, 0.20, 0.22]); leg_ax.set_axis_off()

im = ax.imshow(Hm, aspect='auto', origin='lower', interpolation='nearest',
               vmin=vmin, vmax=vmax, cmap=cmap)

# manual labels axis
ax.set_yticks(np.arange(X.shape[0])); ax.set_yticklabels([])
ax.yaxis.tick_right(); ax.yaxis.set_label_position('right')
ax.set_xticks([]); ax.set_xlabel(''); ax.set_ylabel('')

rax = ax.twinx(); rax.set_ylim(ax.get_ylim()); rax.set_yticks([]); rax.set_xticks([]); rax.set_frame_on(False)
fig.canvas.draw(); renderer = fig.canvas.get_renderer()
x0 = 1.005; tr = blend(rax.transAxes, rax.transData)

def draw_piece(text, x, y, color, weight='normal', style='normal', size=None):
    t = rax.text(x, y, text, transform=tr, ha='left', va='center',
                 fontsize=size or fs_genus, color=color, fontweight=weight, fontstyle=style)
    bb = t.get_window_extent(renderer=renderer)
    dx = bb.width / ax.bbox.width
    return x + dx

for i, raw in enumerate(pretty_labels):
    prefix, core, suffix, mode = split_genus_parenthetical(raw)
    key = raw.strip().lower()
    is_bold = (key in BOLD_SET)
    color = 'black' if is_bold else GRAY_COLOR
    w_all = 'bold' if is_bold else 'normal'
    x = x0
    if mode == 'plain':
        x = draw_piece(raw,   x, i, color, weight=w_all, style='normal')
    elif mode == 'genus':
        if prefix: x = draw_piece(prefix, x, i, color, weight=w_all, style='normal')
        w_gen = 'bold' if is_bold else 'normal'
        x = draw_piece(core,  x, i, color, weight=w_gen, style='italic')
        if suffix: x = draw_piece(suffix, x, i, color, weight=w_all, style='normal')
    elif mode == 'cand':
        x = draw_piece(prefix, x, i, color, weight=w_all, style='italic')
        w_gen = 'bold' if is_bold else 'normal'
        x = draw_piece(core,   x, i, color, weight=w_gen, style='normal')
        if suffix: x = draw_piece(suffix, x, i, color, weight=w_all, style='normal')
    elif mode == 'incertae':
        x = draw_piece(prefix, x, i, color, weight=w_all, style='italic')
        if suffix: x = draw_piece(suffix, x, i, color, weight=w_all, style='normal')

# spacers & site boundaries
for i,sid in enumerate(ordered_cols):
    if sid is None: ax.axvspan(i-0.5, i+0.5, color='lightgray', alpha=0.24, linewidth=0)
for (start,end,bs) in block_bounds[:-1]:
    ax.axvline(end+0.5 + SPACER_COLS_BODYSITE/2.0, color='k', lw=0.7, alpha=0.26)

# top strip + site titles
strip_ax = fig.add_axes([0.19, 0.92, 0.78, 0.020])
strip_ax.imshow(fac_strip, aspect='auto', interpolation='nearest'); strip_ax.set_axis_off()
for (start,end,bs) in block_bounds:
    strip_ax.text((start+end)/2, -1.35, bs, ha='center', va='top',
                  fontsize=fs_site, fontweight=fw_site)

# colorbar
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_ticks(np.linspace(vmin, vmax, 5)); cb.ax.tick_params(labelsize=9)
cb.set_label(f'log10(Relative abundance + {PSEUDO_LOG:.0e})', fontsize=10)

# facility legend
leg_ax.text(0.00, 1.02, 'Facility', ha='left', va='bottom', fontsize=12, fontweight='bold',
            transform=leg_ax.transAxes)
yy = 0.90
for k in facility_order:
    leg_ax.add_patch(Rectangle((0.00, yy-0.08), 0.12, 0.12, color=fac_colors[k], transform=leg_ax.transAxes))
    leg_ax.text(0.15, yy-0.02, k, va='center', ha='left', fontsize=10, transform=leg_ax.transAxes)
    yy -= 0.16
leg_ax.set_xlim(0,1); leg_ax.set_ylim(0,1)

# Save (v3p filenames)
for ext in ('png','pdf','svg'):
    plt.savefig(os.path.join(OUTD, f'genus_heatmap_onepanel_v3p.{ext}'),
                dpi=300 if ext=='png' else None, bbox_inches='tight')

# ================= Font report =================
from matplotlib.text import Text
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
        sample = (t.get_text() or '').replace('\n',' ')[:60]
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

report_path = os.path.join(OUTD, "font_report_genus_heatmap_v3p.txt")
export_font_report(plt.gcf(), report_path, {"output_png": os.path.join(OUTD, "genus_heatmap_onepanel_v3p.png")})
print("Wrote font report:", report_path)
print("Done: genus_heatmap_onepanel_v3p.[png|pdf|svg]")
PY
