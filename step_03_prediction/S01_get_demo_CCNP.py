# -*- coding: utf-8 -*-

SUBLIST_TXT       = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/results/subid_final.txt"                 
DEMO_CSV          = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/code/4th_prediction/demo/CCNP_demo_withoutNAN.csv"    
FD_SUMMARY_CSV    = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/results/head_motion/fd_summary.csv"           
OUT_CSV           = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/code/4th_prediction/demo/demo_age_sex_197.csv"                   


import pandas as pd
import numpy as np
import re
from pathlib import Path

def normalize_colname(c: str) -> str:
    s = str(c).replace("\ufeff", "").strip().lower()
    s = s.replace(" ", "").replace("-", "_")
    return s

def norm_subid(x: str) -> str:
    s = str(x).strip()
    return s if s.startswith("sub-") else f"sub-{s}"

def read_table_smart(path: str) -> pd.DataFrame:
    enc = "utf-8-sig"
    try:
        df = pd.read_csv(path, engine="python", sep=None, encoding=enc)
        if df.shape[1] > 1:
            return df
    except Exception:
        pass
    try:
        df = pd.read_csv(path, engine="python", sep=",", encoding=enc)
        if df.shape[1] > 1:
            return df
    except Exception:
        pass
    try:
        df = pd.read_csv(path, engine="python", sep="\t", encoding=enc)
        if df.shape[1] > 1:
            return df
    except Exception:
        pass
    df = pd.read_csv(path, engine="python", sep=";", encoding=enc)
    return df

def find_subid_column(cols_norm):
    candidates = ["subid_bids", "subid", "subidbids", "sub_id_bids", "sub_id", "participant_id", "subject", "subjectid"]
    for cand in candidates:
        if cand in cols_norm:
            return cand
    for c in cols_norm:
        if ("sub" in c and "bids" in c) or re.fullmatch(r"sub[_]?id", c):
            return c
    return None

def pick_first_existing(cols_norm, *names):
    for n in names:
        if n in cols_norm:
            return n
    return None


def main():

    with open(SUBLIST_TXT, "r", encoding="utf-8") as f:
        subids_raw = [line.strip() for line in f if line.strip()]
    subids = [norm_subid(x) for x in subids_raw]
    order_df = pd.DataFrame({"subid": subids, "_order": range(len(subids))})


    demo = read_table_smart(DEMO_CSV)
    orig_demo_cols = list(demo.columns)
    demo.columns = [normalize_colname(c) for c in demo.columns]

    subcol_demo = find_subid_column(demo.columns)
    if subcol_demo is None:
        raise KeyError(
            f"Current column name (after normalization):{demo.columns.tolist()}; Original column names:{orig_demo_cols}"
        )

    age_col = pick_first_existing(demo.columns, "age")
    sex_col = pick_first_existing(demo.columns, "sex")
    if age_col is None or sex_col is None:
        raise KeyError(
            f"Current column name (after normalization):{demo.columns.tolist()}; Original column names:{orig_demo_cols}"
        )

    demo["subid"] = demo[subcol_demo].astype(str).map(norm_subid)
    demo_small = demo[["subid", age_col, sex_col]].copy()
    demo_small.rename(columns={age_col: "age", sex_col: "sex"}, inplace=True)
    demo_small["age"] = pd.to_numeric(demo_small["age"], errors="coerce")
    demo_small["sex"] = pd.to_numeric(demo_small["sex"], errors="coerce")  # sex

    fd = read_table_smart(FD_SUMMARY_CSV)
    orig_fd_cols = list(fd.columns)
    fd.columns = [normalize_colname(c) for c in fd.columns]

    subcol_fd = find_subid_column(fd.columns)
    if subcol_fd is None:
        if "subid" in fd.columns:
            subcol_fd = "subid"
        else:
            raise KeyError(
                f"Current column name (after normalization):{fd.columns.tolist()}; Original column names:{orig_fd_cols}"
            )

    def find_fd_col(prefix):
        names = [
            f"{prefix}_meanfd", f"{prefix}meanfd",
            f"{prefix}_mean_fd", f"{prefix}mean_fd",
            f"{prefix}_fdmean", f"{prefix}fdmean",
        ]
        for n in names:
            if n in fd.columns:
                return n
        for c in fd.columns:
            if prefix in c and "meanfd" in c:
                return c
        return None

    rest1_col = find_fd_col("rest1")
    rest2_col = find_fd_col("rest2")
    if rest1_col is None and rest2_col is None:
        raise KeyError(
            "The rest1_meanFD / rest2_meanFD (or its variants) columns were not found in the average fd_summary.csv file;"
            f"Current column name (after normalization):{fd.columns.tolist()}; Original column names:{orig_fd_cols}"
        )

    fd["subid"] = fd[subcol_fd].astype(str).map(norm_subid)
    if rest1_col is not None:
        fd[rest1_col] = pd.to_numeric(fd[rest1_col], errors="coerce")
    if rest2_col is not None:
        fd[rest2_col] = pd.to_numeric(fd[rest2_col], errors="coerce")
    fd["meanFD"] = fd[[c for c in [rest1_col, rest2_col] if c is not None]].mean(axis=1, skipna=True)
    fd_small = fd[["subid", "meanFD"]].copy()

    # 3) Merge and output in sublist order
    out = (
        order_df
        .merge(demo_small, on="subid", how="left")
        .merge(fd_small, on="subid", how="left")
        .sort_values("_order")
        .drop(columns=["_order"])
    )
    out = out[["subid", "age", "sex", "meanFD"]]

    n = len(out)
    print(f"[complete] total {n}  subid")
    print(f"missing -> age: {int(out['age'].isna().sum())}, "
          f"sex: {int(out['sex'].isna().sum())}, "
          f"meanFD: {int(out['meanFD'].isna().sum())}")

    uniq_sex = sorted(set(out["sex"].dropna()))
    bad_sex = [v for v in uniq_sex if v not in (1, 2)]
    if bad_sex:
        print(f" sex exists a unique value that is not 1/2: {bad_sex}")

    miss_in_demo = sorted(set(subids) - set(demo_small["subid"]))
    miss_in_fd   = sorted(set(subids) - set(fd_small["subid"]))
    if miss_in_demo:
        print(f"demo Table (Age/Gender) Missing {len(miss_in_demo)}, example {miss_in_demo[:5]}")
    if miss_in_fd:
        print(f"FD table(meanFD) missing {len(miss_in_fd)}, example {miss_in_fd[:5]}")

    out_path = Path(OUT_CSV)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"[save] {out_path.resolve()}")

if __name__ == "__main__":
    main()
