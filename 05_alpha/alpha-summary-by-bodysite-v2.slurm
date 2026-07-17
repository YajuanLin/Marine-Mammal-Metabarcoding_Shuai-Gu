#!/bin/bash
#SBATCH --job-name=alpha-summary-by-bodysite-v2
#SBATCH --output=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-summary-by-bodysite-v2.out
#SBATCH --error=/home/sgu1/sg/formal/merge4tubes/alpha/more/alpha-summary-by-bodysite-v2.err
#SBATCH --time=02:10:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

set -euo pipefail
module load qiime2/2024.2   # provides python + pandas (+scipy present in this env)

BASE=~/sg/formal/merge4tubes
OUT=$BASE/alpha/more
DEPTH=10000  # must match your alpha computation depth
META=$BASE/merged_dolphin_4.txt
TT_TSV=$OUT/Tt/alpha_indices_${DEPTH}_Tt.tsv
SB_TSV=$OUT/Sb/alpha_indices_${DEPTH}_Sb.tsv

SUMMARY_DIR=$OUT/summary
mkdir -p "$SUMMARY_DIR"

python <<'PY'
import os, sys, math, itertools
import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu

BASE=os.path.expanduser("~/sg/formal/merge4tubes")
OUT = os.path.join(BASE, "alpha/more")
SUMMARY_DIR = os.path.join(OUT, "summary")
os.makedirs(SUMMARY_DIR, exist_ok=True)

DEPTH=10000
META = os.path.join(BASE, "merged_dolphin_4.txt")
TT   = os.path.join(OUT, "Tt", f"alpha_indices_{DEPTH}_Tt.tsv")
SB   = os.path.join(OUT, "Sb", f"alpha_indices_{DEPTH}_Sb.tsv")

# --- load metadata (tab-delimited), keep key columns
meta = pd.read_csv(META, sep="\t", dtype=str).fillna("")
need_cols = ["SampleID","BodySite","Species","DolphinID"]
for c in need_cols:
    if c not in meta.columns:
        raise SystemExit(f"Metadata missing column: {c}")
# filter out DolphinID=57 just in case (should already be excluded upstream)
meta = meta[meta["DolphinID"]!="57"]

# --- load alpha tables and tag species (defensive: rely on metadata too)
def load_alpha(path, species_label):
    df = pd.read_csv(path, sep="\t")
    if "SampleID" not in df.columns:
        # accept common alternatives
        for alt in ["sample-id","sampleid","id","sample"]:
            if alt in df.columns:
                df = df.rename(columns={alt:"SampleID"})
                break
        if "SampleID" not in df.columns:
            raise SystemExit(f"Alpha table at {path} needs a SampleID-like column.")
    df["Species_from_table"] = species_label
    return df

frames = []
if os.path.exists(TT): frames.append(load_alpha(TT, "Tt"))
if os.path.exists(SB): frames.append(load_alpha(SB, "Sb"))
if not frames:
    raise SystemExit("No alpha TSVs found.")
alpha = pd.concat(frames, ignore_index=True)

# join to metadata to get BodySite and canonical Species
df = alpha.merge(meta[["SampleID","BodySite","Species"]], on="SampleID", how="left")
# fallback if Species missing in metadata
if "Species_from_table" in df.columns:
    df["Species"] = df["Species"].fillna(df["Species_from_table"])
missing_bs = df["BodySite"].isna().sum()
if missing_bs:
    print(f"WARNING: {missing_bs} samples missing BodySite in metadata; they will be dropped.", file=sys.stderr)
df = df.dropna(subset=["BodySite"])

# metrics present
metrics = ["observed_features","shannon","simpson","pielou_e","chao1"]
for m in metrics:
    if m not in df.columns:
        raise SystemExit(f"Alpha table missing metric column: {m}")
    df[m] = pd.to_numeric(df[m], errors="coerce")

# helper: compute boxplot stats for one series
def box_stats(s: pd.Series) -> pd.Series:
    s = s.dropna().astype(float)
    if s.empty:
        return pd.Series(dict(n=0, mean=float("nan"), median=float("nan"),
                              std=float("nan"), q1=float("nan"), q3=float("nan"),
                              iqr=float("nan"), min=float("nan"), max=float("nan")))
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    return pd.Series(dict(
        n=int(s.shape[0]),
        mean=s.mean(),
        median=s.median(),
        std=s.std(ddof=1) if s.shape[0] > 1 else 0.0,
        q1=q1,
        q3=q3,
        iqr=q3-q1,
        min=s.min(),
        max=s.max(),
    ))

# -------- (A) Summary by BodySite (species combined) --------
rows = []
for m in metrics:
    g = df.groupby("BodySite")[m].apply(box_stats).unstack()
    g.insert(0, "metric", m)
    rows.append(g.reset_index())
summary_combined = pd.concat(rows, ignore_index=True)
summary_combined = summary_combined[["metric","BodySite","n","mean","median","std","q1","q3","iqr","min","max"]]

# -------- (B) Summary by BodySite + Species (Tt and Sb split) --------
rows = []
for m in metrics:
    g = df.groupby(["BodySite","Species"])[m].apply(box_stats).unstack()
    g.insert(0, "metric", m)
    rows.append(g.reset_index())
summary_split_long = pd.concat(rows, ignore_index=True)
summary_split_long = summary_split_long[["metric","BodySite","Species","n","mean","median","std","q1","q3","iqr","min","max"]]

# also pivot to wide so each BodySite shows Tt_* and Sb_* side-by-side per metric
def to_wide(df_long):
    stat_cols = ["n","mean","median","std","q1","q3","iqr","min","max"]
    pieces = []
    for m in metrics:
        sub = df_long[df_long["metric"]==m].copy()
        w = sub.pivot(index="BodySite", columns="Species", values=stat_cols)
        # flatten columns: ('mean','Tt') -> 'Tt_mean'
        w.columns = [f"{sp}_{stat}" for stat,sp in w.columns]
        w.insert(0, "metric", m)
        w = w.reset_index()
        pieces.append(w)
    return pd.concat(pieces, ignore_index=True)

summary_split_wide = to_wide(summary_split_long)

# ---------------- write existing outputs ----------------
summary_combined.to_csv(os.path.join(SUMMARY_DIR, "alpha_summary_by_BodySite_combined.tsv"), sep="\t", index=False)
summary_split_long.to_csv(os.path.join(SUMMARY_DIR, "alpha_summary_by_BodySite_and_Species_long.tsv"), sep="\t", index=False)
summary_split_wide.to_csv(os.path.join(SUMMARY_DIR, "alpha_summary_by_BodySite_and_Species_wide.tsv"), sep="\t", index=False)

# ===================== NEW SECTION: statistical comparisons =====================
def bh_fdr(pvals):
    """Benjamini–Hochberg FDR with numpy only."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = np.arange(1, n+1)
    adj = np.empty(n, dtype=float)
    adj[order] = np.minimum.accumulate((p[order][::-1] * n / ranked[::-1]))[::-1]
    return np.clip(adj, 0, 1)

# (C1) Within-species across BodySites
within_rows = []
for sp, sub_sp in df.groupby("Species"):
    for m in metrics:
        # gather groups by bodysite (only groups with >=1 value)
        groups, labels = [], []
        for bs, sub_bs in sub_sp.groupby("BodySite"):
            vals = sub_bs[m].dropna().values
            if len(vals) > 0:
                groups.append(vals); labels.append(bs)
        # overall Kruskal–Wallis if >=3 groups
        if len(groups) >= 3:
            H, p_kw = kruskal(*groups)
            within_rows.append({
                "Scope":"within_species",
                "Species": sp,
                "BodySite": "ALL",
                "Metric": m,
                "Test": "Kruskal-Wallis",
                "Group1": None, "Group2": None,
                "n_groups": len(groups),
                "stat": H,
                "p_value_raw": p_kw,
                "p_value_adj_BH": p_kw
            })
        # pairwise Mann–Whitney U (two-sided) with BH per (species, metric)
        pairs = list(itertools.combinations(range(len(groups)), 2))
        if pairs:
            pvals, tmp = [], []
            for i,j in pairs:
                g1, g2 = groups[i], groups[j]
                if len(g1) and len(g2):
                    U, p = mannwhitneyu(g1, g2, alternative="two-sided")
                    tmp.append((labels[i], labels[j], len(g1), len(g2), U, p))
                    pvals.append(p)
            if pvals:
                p_adj = bh_fdr(pvals)
                for (bs1, bs2, n1, n2, U, p), pa in zip(tmp, p_adj):
                    within_rows.append({
                        "Scope":"within_species",
                        "Species": sp,
                        "BodySite": f"{bs1} vs {bs2}",
                        "Metric": m,
                        "Test": "Mann-Whitney U",
                        "Group1": bs1, "Group2": bs2,
                        "n1": n1, "n2": n2,
                        "stat": U,
                        "p_value_raw": p,
                        "p_value_adj_BH": pa
                    })

within_tbl = pd.DataFrame(within_rows)
within_path = os.path.join(SUMMARY_DIR, "alpha_tests_within_species_across_bodysites.tsv")
if not within_tbl.empty:
    within_tbl.to_csv(within_path, sep="\t", index=False)

# (C2) Across species (Tt vs Sb) within each BodySite
across_rows = []
species_present = sorted(df["Species"].dropna().unique())
has_tt_sb = set(species_present) >= set(["Tt","Sb"])
for bs, sub_bs in df.groupby("BodySite"):
    if not has_tt_sb:
        continue
    for m in metrics:
        g1 = sub_bs.loc[sub_bs["Species"]=="Tt", m].dropna().values
        g2 = sub_bs.loc[sub_bs["Species"]=="Sb", m].dropna().values
        if len(g1) and len(g2):
            U, p = mannwhitneyu(g1, g2, alternative="two-sided")
            across_rows.append({
                "Scope":"across_species",
                "BodySite": bs,
                "Metric": m,
                "Test": "Mann-Whitney U",
                "Species1": "Tt", "Species2": "Sb",
                "n1": len(g1), "n2": len(g2),
                "stat": U,
                "p_value_raw": p
            })

across_tbl = pd.DataFrame(across_rows)
across_path = os.path.join(SUMMARY_DIR, "alpha_tests_across_species_within_bodysite.tsv")
if not across_tbl.empty:
    # BH per metric across body sites
    out_blocks = []
    for m, sub in across_tbl.groupby("Metric", group_keys=False):
        sub = sub.copy()
        sub["p_value_adj_BH"] = bh_fdr(sub["p_value_raw"].values)
        out_blocks.append(sub)
    across_tbl = pd.concat(out_blocks, ignore_index=True)
    across_tbl.to_csv(across_path, sep="\t", index=False)

# ---------------- final console prints ----------------
print("Wrote:")
for f in ["alpha_summary_by_BodySite_combined.tsv",
          "alpha_summary_by_BodySite_and_Species_long.tsv",
          "alpha_summary_by_BodySite_and_Species_wide.tsv"]:
    print("  -", os.path.join(SUMMARY_DIR, f))
if os.path.exists(within_path):
    print("  -", within_path)
else:
    print("  - (no within-species tests written; insufficient groups)")
if os.path.exists(across_path):
    print("  -", across_path)
else:
    print("  - (no across-species tests written; need both Tt and Sb)")
PY

echo "Done."
