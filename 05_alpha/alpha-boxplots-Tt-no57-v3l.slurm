#!/bin/bash
#SBATCH --job-name=alpha-boxplots-Tt-no57-v3l
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-boxplots-Tt-no57-v3l.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-boxplots-Tt-no57-v3l.err
#SBATCH --time=04:20:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shuai.gu@tamucc.edu

set -euo pipefail
module load qiime2/2024.2   # python + matplotlib + pandas

# --- Paths ---
BASE=~/sg/formal/merge4tubes
OUT=$BASE/alpha/more
META=$BASE/merged_dolphin_4.txt
DEPTH=10000
TT_TSV=$OUT/Tt/alpha_indices_${DEPTH}_Tt.tsv
PLOT_DIR=$OUT/plots
mkdir -p "$PLOT_DIR"

# Clear per-user Matplotlib cache to avoid stale font fallback
rm -rf ~/.cache/matplotlib || true

python <<'PY'
import os, random, sys
from collections import OrderedDict, Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

# ---------- 1) Register user-provided Arial fonts ----------
from matplotlib import font_manager
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

# Rebuild font cache so newly added fonts are discoverable
try:
    font_manager._rebuild()
except Exception:
    pass

# ---------- 2) Force Arial everywhere (axis/ticks/legend/titles/annotations/math) ----------
matplotlib.rcParams.update({
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

import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.font_manager import FontProperties, findfont

# --- Font audit helper ---
def export_font_report(fig, report_path, extra=None):
    lines = []
    lines.append(f"Matplotlib version: {matplotlib.__version__}")
    lines.append(f"rcParams[font.family]: {matplotlib.rcParams.get('font.family')}")
    lines.append(f"rcParams[font.sans-serif]: {matplotlib.rcParams.get('font.sans-serif')}")
    lines.append(f"User-registered Arial font files: {len(FOUND)}")
    for p in FOUND:
        lines.append(f"  - {p}")
    resolved = findfont(FontProperties(family=['Arial','sans-serif']))
    lines.append(f"Resolved font path for 'Arial' request: {resolved}")
    lines.append("")

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

# ---- paths ----
BASE = os.path.expanduser("~/sg/formal/merge4tubes")
OUT  = os.path.join(BASE, "alpha/more")
PLOT_DIR = os.path.join(OUT, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)
META = os.path.join(BASE, "merged_dolphin_4.txt")
DEPTH = 10000
TT_TSV = os.path.join(OUT, "Tt", f"alpha_indices_{DEPTH}_Tt.tsv")

# ---- palettes + normalization ----
sys.path.append('/home/sgu1/sg/formal/merge4tubes/palettes')
from palettes_pub import (BODYSITE_ORDER, normalize_bodysite_series)

# Okabe–Ito palette (no black); VP are close blues; others spaced
BODYSITE_COLOR = {
    "Vaginal":      "#56B4E9",
    "Penile":       "#0072B2",
    "Mammary":      "#009E73",
    "Gastric":      "#E69F00",
    "Rectal":       "#CC79A7",
    "Respiratory":  "#F0E442",
}

metrics = OrderedDict([
    ("chao1",             "Chao1 richness"),
    ("shannon",           "Shannon diversity"),
    ("pielou_e",          "Pielou evenness"),
    ("observed_features", "Observed features"),
    ("simpson",           "Simpson diversity"),
])

# ---- metadata ----
meta = pd.read_csv(META, sep="\t", dtype=str).fillna("")
def pick(colnames, *alts):
    low = {c.lower().replace(" ","").replace("_",""): c for c in colnames}
    for a in alts:
        k = a.lower().replace(" ","").replace("_","")
        if k in low: return low[k]
    return None

sid_col = meta.columns[0]
code_col = pick(meta.columns, "code")
bs_col   = pick(meta.columns, "bodysite","body_site")
sp_col   = pick(meta.columns, "species")
did_col  = pick(meta.columns, "dolphinid","animalid","id")

meta = meta.rename(columns={sid_col:"SampleID", code_col:"Code",
                            bs_col:"BodySite", sp_col:"Species", did_col:"DolphinID"})
for c in ["SampleID","Code","BodySite","Species","DolphinID"]:
    meta[c] = meta[c].astype(str).str.strip()

# Exclude DolphinID=57; keep only Species=Tt
meta = meta[(meta["DolphinID"]!="57") &
            (meta["Species"].str.upper().str.startswith("TT"))].copy()

# Normalize BodySite labels and enforce order
meta["BodySite_pub"] = normalize_bodysite_series(meta["BodySite"], meta["DolphinID"])

# ---- alpha (Tt) ----
alpha = pd.read_csv(TT_TSV, sep="\t")
for m in metrics.keys():
    if m not in alpha.columns:
        raise SystemExit(f"Missing metric column: {m}")

# join + keep canonical sites
df = alpha.merge(meta[["SampleID","BodySite_pub"]], on="SampleID", how="inner")
keep_sites = list(BODYSITE_ORDER)
df = df[df["BodySite_pub"].isin(keep_sites)].copy()
ordered_sites = [s for s in keep_sites if s in df["BodySite_pub"].unique()]
rng = random.Random(42)

def draw_panel(ax, df, metric_key, metric_label):
    data, counts, means, colors = [], [], [], []
    for bs in ordered_sites:
        sub = df[df["BodySite_pub"]==bs]
        vals = pd.to_numeric(sub[metric_key], errors="coerce").dropna().values
        data.append(vals)
        counts.append(len(vals))
        means.append(np.mean(vals) if len(vals) else np.nan)
        colors.append(BODYSITE_COLOR.get(bs, "#BBBBBB"))

    positions = list(range(1, len(ordered_sites)+1))

    bp = ax.boxplot(
        data, positions=positions, widths=0.5,
        patch_artist=True, manage_ticks=False, showfliers=False
    )
    for box, c in zip(bp['boxes'], colors):
        box.set(facecolor=c, edgecolor='black', linewidth=1.2, alpha=0.9)
    for key in ('whiskers','caps'):
        for artist in bp[key]:
            artist.set(color='black', linewidth=1.0)
    for median in bp['medians']:
        median.set(color='black', linewidth=1.6)

    for pos, vals, c in zip(positions, data, colors):
        xs = [pos + rng.uniform(-0.12, 0.12) for _ in vals]
        ax.scatter(xs, vals, s=18, alpha=0.6, c=c,
                   edgecolors='black', linewidths=0.4)

    for pos, mu in zip(positions, means):
        if not np.isnan(mu):
            ax.scatter([pos],[mu], marker='D', s=34, c='black',
                       edgecolors='black', linewidths=0.6, zorder=3)

    xticks = [f"{bs}\n n={n}" for bs, n in zip(ordered_sites, counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xticks, rotation=20)

    ax.set_ylabel(metric_label)
    ax.set_xlabel("Body site")
    ax.grid(False)
    ax.margins(x=0.02)

# ---- make multipanel figure (A–C) and write a font report ----
fig, axes = plt.subplots(1, 3, figsize=(18,6), sharex=False, sharey=False)
triplet = [("chao1","Chao1 richness"),
           ("shannon","Shannon diversity"),
           ("pielou_e","Pielou evenness")]
for ax, (k, label) in zip(axes, triplet):
    draw_panel(ax, df, k, label)

for ax, lab in zip(axes, ["A", "B", "C"]):
    ax.text(-0.05, 1.04, lab, transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="top", ha="left")

fig.tight_layout()
out_png = os.path.join(PLOT_DIR, "alpha_box_multi_Tt_no57_by_bodysite_pub_v3l.png")
fig.savefig(out_png, dpi=300)

report_path = os.path.join(PLOT_DIR, "font_report_v3l.txt")
export_font_report(fig, report_path, {"output_png": out_png})
print("Wrote:", out_png)
print("Wrote font report:", report_path)
plt.close(fig)
PY

echo "Done."
