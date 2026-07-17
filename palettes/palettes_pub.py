# palettes_pub.py — publication-locked palettes + normalization helpers
# Location: ~/sg/formal/merge4tubes/palettes/palettes_pub.py
#
# Global label rules across ALL figures:
# - Blowhole -> Respiratory (EXCEPT DolphinID==57; keep "Blowhole")
# - Gastric stays "Gastric"; ONLY "gastric (mouth)" (and close variants) -> "Oral"
# - Internal placenta -> Fetal wall
# - External placenta -> Maternal wall
# - genital_P -> Penile; genital_V -> Vaginal
# Legend/facet base order: ["Vaginal","Penile","Mammary","Gastric","Rectal","Respiratory"]

# ===============================
# Facility/Code palette (no HI3/HI4)
# ===============================
CODE_COLOR = {
    "HI1": "#1B9E77",
    "HI2": "#D95F02",
    "FL1": "#66A61E",
    "FL2": "#E6AB02",
    "FL3": "#A6761D",
    "FL4": "#1F78B4"
}
CODE_ORDER = ["HI1", "HI2", "FL1", "FL2", "FL3", "FL4"]

# ==============================================
# BodySite palette (final publication labels)
# ==============================================
BODYSITE_COLOR = {
    "Vaginal":      "#3366CC",
    "Penile":       "#6699FF",
    "Mammary":      "#DC3912",
    "Gastric":      "#109618",
    "Rectal":       "#990099",
    "Respiratory":  "#FF9900",  # default rename of Blowhole (except DolphinID==57)
    # Additional sites that may appear in some plots:
    "Oral":         "#2CA02C",  # only for "gastric (mouth)" inputs
    "Fetal wall":   "#8C564B",
    "Maternal wall":"#C49C94",
    "Milk":         "#17BECF"
}

# Publication base order
BODYSITE_ORDER = ["Vaginal", "Penile", "Mammary", "Gastric", "Rectal", "Respiratory"]

# ==============================================
# Normalization / aliasing for messy inputs
# ==============================================
_ALIASES = {
    # Genital variants
    "genital_v": "Vaginal",
    "genital-v": "Vaginal",
    "genital v": "Vaginal",
    "vaginal":   "Vaginal",

    "genital_p": "Penile",
    "genital-p": "Penile",
    "genital p": "Penile",
    "penile":    "Penile",

    # Mammary/Rectal
    "mammary": "Mammary",
    "rectal":  "Rectal",

    # Gastric stays Gastric
    "gastric": "Gastric",

    # ONLY these map to Oral:
    "gastric (mouth)": "Oral",
    "gastric-mouth":   "Oral",
    "gastric_mouth":   "Oral",
    "mouth":           "Oral",
    "oral":            "Oral",

    # Blowhole/Respiratory
    "blowhole":    "Respiratory",   # default; exception handled by DolphinID rule
    "respiratory": "Respiratory",

    # Placenta
    "internal placenta": "Fetal wall",
    "placenta internal": "Fetal wall",
    "placenta_internal": "Fetal wall",
    "external placenta": "Maternal wall",
    "placenta external": "Maternal wall",
    "placenta_external": "Maternal wall",

    # Milk
    "milk": "Milk"
}

def _normalize_single(bodysite_raw: str) -> str:
    if bodysite_raw is None:
        return ""
    key = str(bodysite_raw).strip()
    if key == "":
        return ""
    k = key.lower()
    if k in _ALIASES:
        return _ALIASES[k]
    return key.strip().title()

def normalize_bodysite_with_id(bodysite_raw, dolphin_id) -> str:
    """
    Normalize BodySite with DolphinID exception:
    - If bodysite is 'blowhole' and dolphin_id == 57 -> keep 'Blowhole'
    - Else apply standard mapping (Blowhole->Respiratory, gastric(mouth)->Oral, etc.)
    """
    if bodysite_raw is not None and str(bodysite_raw).strip().lower() == "blowhole":
        try:
            if int(dolphin_id) == 57:
                return "Blowhole"
        except Exception:
            pass
    return _normalize_single(bodysite_raw)

def normalize_bodysite_series(bodysite_series, dolphin_id_series):
    """Vectorized convenience for pandas."""
    import pandas as pd
    return pd.Series(
        (normalize_bodysite_with_id(b, d) for b, d in zip(bodysite_series, dolphin_id_series)),
        index=bodysite_series.index
    )

# ===============
# Legend helpers
# ===============
def present_legend_handles(ax, base_order, palette_dict, title=None, frameon=False, loc="best"):
    """
    Build legend handles using 'base_order', but include only labels that are present.
    If present labels not in base_order (e.g., Fetal/Maternal wall, Milk, Oral) exist,
    they are appended after base_order in sorted order so they still appear.
    """
    import matplotlib.patches as mpatches
    present = set()
    for artist in ax.get_children():
        if hasattr(artist, 'get_label'):
            lab = artist.get_label()
            if lab and not lab.startswith('_'):
                present.add(lab)
    ordered = [lab for lab in base_order if lab in present]
    extras = sorted(present.difference(base_order))
    ordered.extend(extras)
    handles = [mpatches.Patch(color=palette_dict.get(l, "#BBBBBB"), label=l) for l in ordered]
    if handles:
        return ax.legend(handles=handles, title=title, frameon=frameon, loc=loc)
    return None

def map_colors(keys, palette_dict, fallback="#BBBBBB"):
    """Return a list of hex colors for the given keys based on the provided palette."""
    return [palette_dict.get(k, fallback) for k in keys]
