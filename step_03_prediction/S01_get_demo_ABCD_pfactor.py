# -*- coding: utf-8 -*-

import os
import re
from pathlib import Path
import pandas as pd

proj_dir = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD")
data_dir = proj_dir / "demo"
out_dir  = proj_dir / "code/5th_prediction/prepare_data"
out_dir.mkdir(parents=True, exist_ok=True)  

demo_csv       = data_dir / "abcd_p_demo_baseline.csv"         # sex/race/parent edu/income/family size
long_csv       = data_dir / "abcd_y_lt_baseline.csv"           # age/interview date/site_id/rel_family_id
mri_info_csv   = data_dir / "mri_y_qc_motion_baseline.csv"     # motion
adm_info_csv   = data_dir / "mri_y_adm_info_baseline.csv"      # scanner
pfactor_csv    = data_dir / "Pfactor_score_wx.csv"             # p_factor: General, Ext, ADHD, Int

sublist_txt    = data_dir / "has_smooth_data_new.txt"


motion2run_csv = data_dir / "ABCD_MotionInfo_n5959_corrected_withTPchecked.csv"

out_main_csv   = out_dir / "abcd_final_selected_columns.csv"
out_cov_csv    = out_dir / "abcd_covariates.csv"
out_pf_csv     = out_dir / "abcd_pfactor_table.csv"  

def canonicalize_ndar(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip().str.upper()
    return s.str.replace(r"^NDAR[_-]?INV", "NDAR_INV", regex=True)

def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns=lambda c: re.sub(r"\ufeff", "", str(c)).strip())

def ensure_subject_id(df: pd.DataFrame, name="df") -> pd.DataFrame:
    id_candidates = ["src_subject_id", "subject_id", "subjectkey", "id", "subList", "SUBLIST"]
    hit = next((c for c in id_candidates if c in df.columns), None)
    if hit is None:
        lower_map = {c.lower(): c for c in df.columns}
        for c in [x.lower() for x in id_candidates]:
            if c in lower_map:
                hit = lower_map[c]; break
    if hit is None:
        raise ValueError(f"{name}: Subject ID column not found; acceptable column names: {id_candidates};actual column:{list(df.columns)[:10]}...")
    if hit != "subject_id":
        df = df.rename(columns={hit: "subject_id"})
    df["subject_id"] = canonicalize_ndar(df["subject_id"])
    return df

def read_csv_with_id(path: Path, name: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    df = clean_columns(df)
    return ensure_subject_id(df, name=name)

def to_numeric_if_needed(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.strip().replace({"": pd.NA}), errors="coerce")

def make_site_id_numeric(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    s_num = pd.to_numeric(s, errors="coerce")
    if s_num.notna().all():
        return s_num.astype("Int64")

    extracted = s.str.extract(r"(\d+)")[0]
    s_num2 = pd.to_numeric(extracted, errors="coerce")
    if s_num2.notna().all():
        return s_num2.astype("Int64")

    codes, _ = pd.factorize(s, sort=True)
    return pd.Series(codes + 1, index=series.index, dtype="Int64")


demoT = read_csv_with_id(demo_csv, name="demo")
longT = read_csv_with_id(long_csv, name="longitudinal")
mriT  = read_csv_with_id(mri_info_csv, name="mri_info")
admT  = read_csv_with_id(adm_info_csv, name="adm_info")
pfT   = read_csv_with_id(pfactor_csv, name="p_factor")
m2r   = read_csv_with_id(motion2run_csv, name="motion_2run")


with open(sublist_txt, "r", encoding="utf-8") as f:
    sub_ids = [ln.strip() for ln in f if ln.strip()]
subList = set(canonicalize_ndar(pd.Series(sub_ids)).tolist())


# demo: sex, race
demo_cols_map = {
    "demo_sex_v2": "sex",
    "race_ethnicity": "race",
}
missing_demo = [c for c in demo_cols_map if c not in demoT.columns]
if missing_demo:
    raise KeyError(f"demo missing column: {missing_demo}")
used_demo = demoT[["subject_id"] + list(demo_cols_map.keys())].rename(columns=demo_cols_map)


long_cols_map = {
    "interview_age": "age",           
    "interview_date": "interview_date",
    "site_id_l": "site_id",
    "rel_family_id": "real_family_id",
}
missing_long = [c for c in long_cols_map if c not in longT.columns]
if missing_long:
    raise KeyError(f"longitudinal missing column: {missing_long}")
used_long = longT[["subject_id"] + list(long_cols_map.keys())].rename(columns=long_cols_map)

# mri: mean_motion
mri_cols_map = {"rsfmri_meanmotion": "mean_motion"}
missing_mri = [c for c in mri_cols_map if c not in mriT.columns]
if missing_mri:
    raise KeyError(f"mri_info missing column: {missing_mri}")
used_mri = mriT[["subject_id"] + list(mri_cols_map.keys())].rename(columns=mri_cols_map)

# adm: scanner
adm_cols_map = {"mri_info_manufacturer": "scanner"}
missing_adm = [c for c in adm_cols_map if c not in admT.columns]
if missing_adm:
    raise KeyError(f"adm_info missing column: {missing_adm}")
used_adm = admT[["subject_id"] + list(adm_cols_map.keys())].rename(columns=adm_cols_map)

# p_factor: General/Ext/ADHD/Int
pf_cols_map = {"General": "General", "Ext": "Ext", "ADHD": "ADHD", "Int": "Int"}
missing_pf = [c for c in pf_cols_map if c not in pfT.columns]
if missing_pf:
    raise KeyError(f"p_factor missing column: {missing_pf}")
used_pf = pfT[["subject_id"] + list(pf_cols_map.keys())].rename(columns=pf_cols_map)

# motion2run: meanFD + run1/run2 +  FDratio_lessthan_02
m2r_needed = [
    "meanFD",
    "meanFD_run1", "FD02ratio_run1", "QCvalid_run1",
    "meanFD_run2", "FD02ratio_run2", "QCvalid_run2",
    "FDratio_lessthan_02" 
]
missing_m2r = [c for c in m2r_needed if c not in m2r.columns]
if missing_m2r:
    raise KeyError(f"motion_2run missing column: {missing_m2r}")
used_m2r = m2r[["subject_id"] + m2r_needed].copy()

# ----------------  subList  ----------------
def keep_sublist(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["subject_id"].isin(subList)].copy()

used_demo = keep_sublist(used_demo)
used_long = keep_sublist(used_long)
used_mri  = keep_sublist(used_mri)
used_adm  = keep_sublist(used_adm)
used_pf   = keep_sublist(used_pf)
used_m2r  = keep_sublist(used_m2r)


merged = (
    used_long
    .merge(used_demo, on="subject_id", how="inner")
    .merge(used_mri,  on="subject_id", how="inner")
    .merge(used_adm,  on="subject_id", how="inner")
    .merge(used_pf,   on="subject_id", how="inner")
    .merge(used_m2r,  on="subject_id", how="inner")
)


for c in ["age", "sex", "mean_motion", "meanFD", "FDratio_lessthan_02", "General", "Ext", "ADHD", "Int"]:
    if c in merged.columns:
        merged[c] = to_numeric_if_needed(merged[c])

# delete sex == 3
if "sex" in merged.columns:
    merged = merged.loc[~(merged["sex"] == 3)].copy()

# meanFD <= 0.5 & FDratio_lessthan_02 > 0.4
if not {"meanFD", "FDratio_lessthan_02"}.issubset(merged.columns):
    raise KeyError("missing meanFD or FDratio_lessthan_02")
merged = merged.loc[
    (merged["meanFD"] <= 0.5) &
    (merged["FDratio_lessthan_02"] > 0.4)
].copy()


required_numeric = ["sex", "age", "General", "Ext", "ADHD", "Int", "meanFD"]
merged = merged.dropna(subset=[c for c in required_numeric if c in merged.columns], how="any").copy()


merged = merged.drop_duplicates(subset=["subject_id"]).reset_index(drop=True)


cols_to_drop = [c for c in merged.columns if c.lower() == "site"]
if cols_to_drop:
    merged = merged.drop(columns=cols_to_drop)


final_cols = [
    "subject_id", "age", "interview_date", "site_id", "real_family_id",
    "sex", "race", "mean_motion", "scanner", "meanFD",
    "General", "Ext", "ADHD", "Int",
    "meanFD_run1", "FD02ratio_run1", "QCvalid_run1",
    "meanFD_run2", "FD02ratio_run2", "QCvalid_run2",
    "FDratio_lessthan_02" 
missing_final = [c for c in final_cols if c not in merged.columns]
if missing_final:
    raise KeyError(f"Finally, export the missing columns: {missing_final}")

main_table = merged[final_cols].copy()
main_table.to_csv(out_main_csv, index=False)

cov = main_table[["subject_id", "age", "sex", "meanFD", "site_id"]].copy()
cov["site_id"] = make_site_id_numeric(cov["site_id"])
cov.to_csv(out_cov_csv, index=False)

pf_table = main_table[["subject_id", "age", "interview_date", "site_id", "sex", "meanFD", "General", "Ext", "ADHD", "Int"]].copy()
pf_table.to_csv(out_pf_csv, index=False)

print(f"[SAVE] table: {out_main_csv}")
print(f"[SAVE] Covariate table: {out_cov_csv}")
print(f"[SAVE] p-factor table: {out_pf_csv}")
