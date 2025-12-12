# -*- coding: utf-8 -*-

SUBID_TXT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results/subid_final.txt"
FD_CSV = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results/head_motion/fd_summary.csv"
PARTICIPANTS_TSV = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/demography/study-PNC_desc-participants.tsv"
OUT_CSV = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/code/4th_prediction/demo/demography_age_sex_PNC.csv"
# ======================================

import pandas as pd
import numpy as np
from pathlib import Path

def norm_subid(x: str) -> str:
    s = str(x).strip()
    return s if s.startswith("sub-") else f"sub-{s}"

def map_sex(v):
    if pd.isna(v):
        return np.nan
    s = str(v).strip().lower()
    if s in {"male", "m", "1"}:
        return 1
    if s in {"female", "f", "2"}:
        return 2
    return np.nan

def main():
    subid_path = Path(SUBID_TXT)
    with subid_path.open("r", encoding="utf-8") as f:
        subids = [line.strip() for line in f if line.strip()]
    order_df = pd.DataFrame({"subid": subids, "_order": range(len(subids))})

    fd = pd.read_csv(FD_CSV, engine="python", sep=None)  
    fd.rename(columns={c: c.strip() for c in fd.columns}, inplace=True)
    if "subid" not in fd.columns:
        raise KeyError("fd_summary.csv missing column 'subid'")
    if "rest1_meanFD" not in fd.columns:
        raise KeyError("fd_summary.csv missing column 'rest1_meanFD'")
    fd["subid"] = fd["subid"].astype(str).map(norm_subid)
    fd_small = fd[["subid", "rest1_meanFD"]].copy()
    fd_small.rename(columns={"rest1_meanFD": "meanFD"}, inplace=True)
    fd_small["meanFD"] = pd.to_numeric(fd_small["meanFD"], errors="coerce")

    parts = pd.read_csv(PARTICIPANTS_TSV, sep="\t")
    parts.rename(columns={c: c.strip() for c in parts.columns}, inplace=True)
    for col in ["participant_id", "age", "sex"]:
        if col not in parts.columns:
            raise KeyError(f"study-PNC_desc-participants.tsv missing column '{col}'")
    parts["subid"] = parts["participant_id"].astype(str).map(norm_subid)
    parts_small = parts[["subid", "age", "sex"]].copy()
    parts_small["age"] = pd.to_numeric(parts_small["age"], errors="coerce")
    parts_small["sex"] = parts_small["sex"].apply(map_sex)

    demo = (
        order_df
        .merge(parts_small, on="subid", how="left")
        .merge(fd_small, on="subid", how="left")
        .sort_values("_order")
        .drop(columns=["_order"])
    )
    demo = demo[["subid", "age", "sex", "meanFD"]]

    # 5) missing 
    total = len(demo)
    miss_age = int(demo["age"].isna().sum())
    miss_sex = int(demo["sex"].isna().sum())
    miss_fd  = int(demo["meanFD"].isna().sum())
    print(f"[completed] total {total}  subid; missing -> age:{miss_age}, sex:{miss_sex}, meanFD:{miss_fd}")

    in_subid = set(subids)
    in_parts = set(parts_small["subid"])
    in_fd    = set(fd_small["subid"])
    miss_in_parts = sorted(in_subid - in_parts)
    miss_in_fd    = sorted(in_subid - in_fd)
    if miss_in_parts:
        print(f"participants missing subid {len(miss_in_parts)}, example:{miss_in_parts[:5]}")
    if miss_in_fd:
        print(f"fd_summary missing subid {len(miss_in_fd)}, example: {miss_in_fd[:5]}")

    out_path = Path(OUT_CSV)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    demo.to_csv(out_path, index=False)
    print(f"[save] {out_path.resolve()}")

if __name__ == "__main__":
    main()
