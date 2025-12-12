#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
import numpy as np

FD_SUMMARY = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/head_motion/fd_summary_addframe.csv"
SUBLIST_TXT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/sublist_final_532.txt"
# DEMO_XLSX  = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/demography/EFNY_demo_532.xlsx"
DEMO_XLSX  = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/demography/EFNY_demo_deleteADHD.xlsx"

OUT_NAME   = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/4th_prediction/s01_preparedata/subid_deleteADHD.csv"

# ====== Naming rules for the frame column & Expected number of frames ======
FRAME_SUFFIX = "_nFrames"   # rest1_keep -> rest1_nFrames
FRAME_EXPECT = 179          #  179 
# ====================================================

def _norm_subid(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip()
    return s[4:] if s.lower().startswith("sub-") else s

def _sex_to_code(x):
    if x is None:
        return pd.NA
    s = str(x).strip().lower()
    if s in {"男", "m", "male"}:
        return 1
    if s in {"女", "f", "female"}:
        return 2
    try:
        v = int(float(s))
        if v in (1, 2):
            return v
    except Exception:
        pass
    return pd.NA

def main():
    fd_summary_path = Path(FD_SUMMARY)
    sublist_path    = Path(SUBLIST_TXT)
    demo_path       = Path(DEMO_XLSX)

    # --- sublist ---
    subs = []
    with sublist_path.open("r", encoding="utf-8") as f:
        for line in f:
            t = line.strip()
            if t:
                subs.append(t)
    sublist_norm = {_norm_subid(s) for s in subs}

    # --- FD summary ---
    df_fd = pd.read_csv(fd_summary_path)
    df_fd["subid_norm"] = df_fd["subid"].astype(str).map(_norm_subid)
    df_fd = df_fd[df_fd["subid_norm"].isin(sublist_norm)].copy()

    # --- mean FD: only use keep==1 and frame==FRAME_EXPECT ---
    # Automatically find restX_keep and the corresponding meanFD and frame columns
    keep_cols = [c for c in df_fd.columns if c.startswith("rest") and c.endswith("_keep")]
    mean_cols = [c.replace("_keep", "_meanFD") for c in keep_cols]

    # rest1_keep -> rest1_nFrames
    frame_cols = []
    for k_col in keep_cols:
        if k_col.endswith("_keep"):
            base = k_col[:-5]  # remove "_keep"
        else:
            base = k_col
        cand = base + FRAME_SUFFIX
        frame_cols.append(cand if cand in df_fd.columns else None)

    def row_mean_fd(row):
        vals = []
        for k_col, m_col, f_col in zip(keep_cols, mean_cols, frame_cols):
            
            keep = row.get(k_col, 0)
            if pd.isna(keep):
                keep = 0
            try:
                keep = int(keep)
            except Exception:
                keep = 0

            # meanFD
            m = row.get(m_col, np.nan)

            # frame number
            if f_col is not None and f_col in row.index:
                frames = row.get(f_col, np.nan)
                try:
                    frames = int(frames) if pd.notna(frames) else np.nan
                except Exception:
                    frames = np.nan
            else:
            
                frames = np.nan

            # keep==1 and frame==FRAME_EXPECT and meanFD no NaN
            if (
                keep == 1
                and pd.notna(m)
                and pd.notna(frames)
                and int(frames) == FRAME_EXPECT
            ):
                vals.append(float(m))

        return float(np.mean(vals)) if len(vals) >= 1 else np.nan

    df_fd["mean_FD"] = df_fd.apply(row_mean_fd, axis=1)

    # --- demo (Excel: columns are subid, age, sex) ---
    df_demo = pd.read_excel(demo_path, dtype={"subid": str})
    need_cols = ["subid", "sex", "age"]
    missing = [c for c in need_cols if c not in df_demo.columns]
    if missing:
        raise ValueError(f"demo missing column: {missing}")

    df_demo = df_demo[need_cols].copy()
    # If the same subid has multiple rows, keep only the last row.
    df_demo = df_demo.drop_duplicates(subset=["subid"], keep="last")

    df_demo["subid_norm"] = df_demo["subid"].astype(str).map(_norm_subid)
    df_demo["sex"] = df_demo["sex"].map(_sex_to_code).astype("Int64")
    df_demo["age"] = pd.to_numeric(df_demo["age"], errors="coerce")

    # --- merge ---
    merged = (
        df_fd[["subid", "subid_norm", "mean_FD"]]
        .merge(
            df_demo[["subid_norm", "age", "sex"]],
            on="subid_norm",
            how="left",
            validate="m:1",
        )
    )

    out_df = merged[["subid", "age", "sex", "mean_FD"]].copy()

    # --- write ---
    out_path = Path(OUT_NAME)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[OK] Wrote {out_path} with {len(out_df)} rows.")

if __name__ == "__main__":
    main()