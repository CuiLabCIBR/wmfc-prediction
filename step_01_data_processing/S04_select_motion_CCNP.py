#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
from pathlib import Path
from typing import Tuple, Optional, List, Dict
import re
import sys
import numpy as np
import pandas as pd

SUBLIST = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/batchcode/sublist_all.txt"

# fMRIPrep dir
BASES = [
    "/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/CCNP/fmriprep",
    "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/results/fmriprep",
]

TASK = "rest"
MAX_RUNS = 2
SESSION = "ses-01" 
OUTPUT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/results/head_motion/fd_summary.csv"
PASSED_TXT_OUT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/CCNP/results/subid_final.txt"    


RUN_PATTERN = re.compile(r"run[-_]?(\d+)", re.IGNORECASE)
TASK_PATTERN_FMT = r"task-{}"               

def _normalize_subid(raw: str) -> Tuple[str, str]:
    s = raw.strip()
    if not s:
        return ("", "")
    if s.startswith("sub-"):
        return (s, s)
    return (s, f"sub-{s}")

def _compute_fd_metrics(tsv_path: Path) -> Tuple[Optional[float], Optional[float], int]:
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"[WARN] Failed to read {tsv_path}: {e}")
        return (None, None, 0)

    if "framewise_displacement" not in df.columns:
        print(f"[WARN] No 'framewise_displacement' in {tsv_path}")
        return (None, None, 0)

    fd = pd.to_numeric(df["framewise_displacement"], errors="coerce").dropna()
    n = int(fd.shape[0])
    if n == 0:
        return (None, None, 0)

    mean_fd = float(fd.mean())
    prop_leq = float((fd <= 0.2).mean())
    return (mean_fd, prop_leq, n)

def _find_confounds_for_sub(
    subid: str,
    bases: List[Path],
    task: str,
    max_runs: int,
    session: Optional[str] = None,
) -> Dict[int, Path]:
    task_token = TASK_PATTERN_FMT.format(task.lower())
    as_is, with_sub = _normalize_subid(subid)
    subj_variants = [p for p in [as_is, with_sub] if p]

    found: Dict[int, Path] = {}
    no_run_bucket: List[Path] = []

    for base in bases:
        for subj in subj_variants:
            root = (base / subj)
            if not root.exists():
                continue

            search_root = root
            if session:
                sr = root / session / "func"
                if not sr.exists():
                    continue
                search_root = sr

            for path in search_root.rglob("*desc-confounds_timeseries.tsv"):
            
                if "func" not in [p.name for p in path.parents]:
                    continue

                name = path.name.lower()
                if task_token not in name:
                    continue

                lower_name = path.name.lower()
                if as_is and as_is.lower() not in lower_name and with_sub.lower() not in lower_name:
                    continue

                m = RUN_PATTERN.search(name)
                run_num: Optional[int]
                if m:
                    try:
                        run_num = int(m.group(1))
                    except Exception:
                        run_num = None
                else:
                    run_num = None

                if run_num is None:
                    no_run_bucket.append(path)
                else:
                    if 1 <= run_num <= max_runs and run_num not in found:
                        found[run_num] = path

        if len(found) >= max_runs:
            break

    if no_run_bucket:
        for path in no_run_bucket:
            for r in range(1, max_runs + 1):
                if r not in found:
                    found[r] = path
                    break
            if len(found) >= max_runs:
                break

    return found

def summarize_subject(
    subid: str,
    bases: List[Path],
    task: str,
    max_runs: int,
    session: Optional[str] = None,
):
    runs = _find_confounds_for_sub(subid, bases, task, max_runs, session=session)
    result = {"subid": subid}
    keeps = []
    for r in range(1, max_runs + 1):
        mean_fd = np.nan
        prop_leq = np.nan
        keep = 0
        if r in runs:
            m, p, n = _compute_fd_metrics(runs[r])
            if m is not None and p is not None and n > 0:
                mean_fd = m
                prop_leq = p

                if (m <= 0.5) and (p >= 0.40):
                    keep = 1
        result[f"rest{r}_meanFD"] = mean_fd
        result[f"rest{r}_propFD_leq_0p2"] = prop_leq
        result[f"rest{r}_keep"] = keep
        keeps.append(keep)
    result["check_has2runs"] = 1 if sum(keeps) >= 2 else 0
    return result

def _to_abs(p: Path, base: Path) -> Path:
    return p if p.is_absolute() else (base / p).resolve()

def main():
    script_dir = Path(__file__).resolve().parent

    # sublist
    sublist_path = _to_abs(Path(SUBLIST), script_dir)
    if not sublist_path.exists():
        print(f"[ERR] sublist not found: {sublist_path}", file=sys.stderr)
        sys.exit(1)

    # bases
    bases: List[Path] = []
    for b in BASES:
        p = _to_abs(Path(b), script_dir)
        if not p.exists():
            print(f"[WARN] base path not found: {p}")
        bases.append(p)

    # subjects
    subs: List[str] = []
    with sublist_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            subs.append(line)

    if not subs:
        print("[ERR] empty sublist.", file=sys.stderr)
        sys.exit(1)

    # summarize
    rows = []
    for sub in subs:
        rows.append(summarize_subject(sub, bases, TASK, MAX_RUNS, session=SESSION))

    # dataframe & columns
    df = pd.DataFrame(rows)
    cols = ["subid"]
    for r in range(1, MAX_RUNS + 1):
        cols += [f"rest{r}_meanFD", f"rest{r}_propFD_leq_0p2", f"rest{r}_keep"]
    cols += ["check_has2runs"]
    df = df[cols]


    try:
        df = df.sort_values("subid")
    except Exception:
        pass


    out_csv = _to_abs(Path(OUTPUT), script_dir)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False, float_format="%.6f", encoding="utf-8")
    print(f"[OK] Wrote summary CSV: {out_csv} ({len(df)} rows)")


    keep_cols = [c for c in df.columns if re.match(r"rest\d+_keep$", c)]
    has_r1r2 = ("rest1_keep" in df.columns) and ("rest2_keep" in df.columns)

    if has_r1r2:
        mask = (df["rest1_keep"] == 1) & (df["rest2_keep"] == 1)
        criterion_used = "both rest1_keep==1 and rest2_keep==1"
    else:
        k = (
            df[keep_cols]
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0)
            .ge(1)
            .sum(axis=1)
        )
        mask = k >= 2
        criterion_used = "at least 2 runs keep==1 (fallback)"

    passed_subids = (
        df.loc[mask, "subid"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
    )

    out_txt = _to_abs(Path(PASSED_TXT_OUT), script_dir)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    passed_subids.to_csv(out_txt, index=False, header=False, encoding="utf-8")
    print(f"[OK] Wrote passed subids TXT: {out_txt} (n={len(passed_subids)})")
    print(f"[INFO] Selection criterion: {criterion_used}")

if __name__ == "__main__":
    main()
