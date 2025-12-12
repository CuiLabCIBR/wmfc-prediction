# -*- coding: utf-8 -*-

from pathlib import Path
import re
import pandas as pd

# ========= data path =========
proj_dir = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD")
data_dir = proj_dir / "demo"
out_dir  = proj_dir / "code/5th_prediction/prepare_data"
out_dir.mkdir(parents=True, exist_ok=True)

motion_csv = data_dir / "ABCD_MotionInfo_n5959_corrected_withTPchecked.csv"
merged_csv = data_dir / "abcd_mergedTable_with_motion2run_deNaN_deSex3_n5829_from_n5959.csv"

out_cov_csv = out_dir / "abcd_covariates_meanFD_cognition.csv"
out_nih_csv = out_dir / "abcd_meanFD_cognition.csv"


# ========= functions =========
def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns=lambda c: re.sub(r"\ufeff", "", str(c)).strip())

def canonicalize_ndar(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip().str.upper()
    # Standardized to the NDAR_INV prefix
    return s.str.replace(r"^NDAR[_-]?INV", "NDAR_INV", regex=True)

def ensure_subject_id(df: pd.DataFrame, name="df") -> pd.DataFrame:
    #  ID 
    candidates = ["subject_id", "src_subject_id", "subjectkey", "subid", "id", "subList", "SUBLIST"]
    col = next((c for c in candidates if c in df.columns), None)
    if col is None:
        lower_map = {c.lower(): c for c in df.columns}
        for c in candidates:
            if c.lower() in lower_map:
                col = lower_map[c.lower()]
                break
    if col is None:
        raise KeyError(f"{name}: The subject ID column was not found; it must contain one of the following:{candidates}; actual column{list(df.columns)[:12]}")
    if col != "subject_id":
        df = df.rename(columns={col: "subject_id"})
    df["subject_id"] = canonicalize_ndar(df["subject_id"])
    return df

def read_csv_with_id(path: Path, name: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    df = clean_columns(df)
    df = ensure_subject_id(df, name=name)
    return df

def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.strip().replace({"": pd.NA}), errors="coerce")

def site_digits_only(series: pd.Series) -> pd.Series:
    """
    Keep only the numbers in site_id (e.g., 'site22' -> 22; if there are no numbers -> NaN).
    """
    extracted = series.astype(str).str.extract(r"(\d+)")[0]
    return pd.to_numeric(extracted, errors="coerce").astype("Int64")


# ========= read and screen motion =========
motion = read_csv_with_id(motion_csv, name="motion")
needed_motion_cols = ["meanFD", "FDratio_lessthan_02"]
missing = [c for c in needed_motion_cols if c not in motion.columns]
if missing:
    raise KeyError(f"motion file missing column: {missing}")

motion["meanFD"] = to_num(motion["meanFD"])
motion["FDratio_lessthan_02"] = to_num(motion["FDratio_lessthan_02"])

motion_filtered = motion.loc[
    (motion["meanFD"] <= 0.5) &
    (motion["FDratio_lessthan_02"] >= 0.4)
][["subject_id", "meanFD"]].dropna(subset=["meanFD"]).drop_duplicates("subject_id")

print(f"[INFO] motion number of participants selected: {len(motion_filtered)}")


# ========= Read in the merged main table =========
merged = read_csv_with_id(merged_csv, name="merged")
# Columns involved in the output requirements
nih_cols = [
    "nihtbx_picvocab","nihtbx_pic","nihtbx_flanker","nihtbx_list",
    "nihtbx_cardsort","nihtbx_pattern","nihtbx_reading",
    "nihtbx_fluidcomp","nihtbx_cryst","nihtbx_totalcomp",
]
base_needed = ["subject_id", "age", "site_id", "sex"]
need_all = base_needed + nih_cols
missing_merged = [c for c in need_all if c not in merged.columns]
if missing_merged:
    raise KeyError(f"merged file missing column: {missing_merged}")

# Keep only the columns you need and remove duplicates.
merged_use = merged[["subject_id", "age", "site_id", "sex"] + nih_cols].drop_duplicates("subject_id")


# ========= Merge (only retain subjects who were selected through motion and are present in the merged list) =========
df = merged_use.merge(motion_filtered, on="subject_id", how="inner")
print(f"[INFO] Number of subjects after merging: {len(df)}")

for c in ["age", "sex", "meanFD"]:
    if c in df.columns:
        df[c] = to_num(df[c])

# site 
df["site_id"] = site_digits_only(df["site_id"])

# ========= export cov =========
cov_cols = ["subject_id", "age", "sex", "meanFD", "site_id"]
cov = df[cov_cols].copy()
cov.to_csv(out_cov_csv, index=False)
print(f"[SAVE] Covariate table: {out_cov_csv}")

# ========= export nih =========
nih_out_cols = ["subject_id", "age", "site_id", "sex"] + nih_cols + ["meanFD"]
nih_out = df[nih_out_cols].copy()
nih_out.to_csv(out_nih_csv, index=False)
print(f"[SAVE] NIH toolbox table: {out_nih_csv}")
