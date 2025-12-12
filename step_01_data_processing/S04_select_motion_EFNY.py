#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Set

BASES = [
    "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/fmriprep",
]

TASK = "rest"
MAX_RUNS = 4
OUTPUT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/head_motion/fd_summary_addframe.csv"

RUN_PATTERN = re.compile(r"run[-_]?(\d+)", re.IGNORECASE)
TASK_PATTERN_FMT = r"task-{}"


def _normalize_subid(raw: str) -> Tuple[str, str]:
    s = raw.strip()
    if not s:
        return ("", "")
    if s.startswith("sub-"):
        return (s, s)
    return (s, f"sub-{s}")


def _enumerate_subjects(bases: List[Path]) -> List[str]:
    subs: Set[str] = set()
    for base in bases:
        if not base.exists():
            continue
        for p in base.iterdir():
            if not p.is_dir():
                continue
            name = p.name.strip()
            if not name:
                continue

            if name.startswith("sub-"):
                subs.add(name)

    return sorted(subs)


def _find_confounds_for_sub(
    subid: str,
    bases: List[Path],
    task: str,
    max_runs: int,
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
            for path in root.rglob("*desc-confounds_timeseries.tsv"):

                if "func" not in [p.name for p in path.parents]:
                    continue
                name = path.name.lower()
                if task_token not in name:
                    continue

                lower_name = path.name.lower()
                if as_is and (as_is.lower() not in lower_name) and (with_sub.lower() not in lower_name):
                    continue

                m = RUN_PATTERN.search(name)
                run_num = None
                if m:
                    try:
                        run_num = int(m.group(1))
                    except Exception:
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


def _compute_fd_metrics(tsv_path: Path) -> Tuple[Optional[float], Optional[float], int]:

    # from 'confounds tsv' read framewise_displacement column：

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


def summarize_subject(
    subid: str,
    bases: List[Path],
    task: str,
    max_runs: int,
):
    runs = _find_confounds_for_sub(subid, bases, task, max_runs)
    result = {"subid": subid}
    keeps = []
    kept_fds = []

    for r in range(1, max_runs + 1):
        mean_fd = np.nan
        prop_leq = np.nan
        n_frames = np.nan    # check frame num
        keep = 0
        if r in runs:
            m, p, n = _compute_fd_metrics(runs[r])
            if m is not None and p is not None and n > 0:
                mean_fd = m
                prop_leq = p
                n_frames = n
                # select FD：meanFD ≤ 0.5 and prop(FD≤0.2) ≥ 0.40
                if (m <= 0.5) and (p >= 0.40):
                    keep = 1
                    kept_fds.append(m)
        result[f"rest{r}_meanFD"] = mean_fd
        result[f"rest{r}_nFrames"] = n_frames           
        result[f"rest{r}_propFD_leq_0p2"] = prop_leq
        result[f"rest{r}_keep"] = keep
        keeps.append(keep)

    kept_count = int(sum(keeps))
    result["kept_runs_count"] = kept_count
    result["check_has2runs"] = 1 if kept_count >= 2 else 0
    result["meanFD_kept_runs"] = float(np.mean(kept_fds)) if kept_count >= 2 else np.nan
    return result


def main():
    script_dir = Path(__file__).resolve().parent

    # base path
    bases: List[Path] = []
    for b in BASES:
        p = Path(b)
        if not p.is_absolute():
            p = (script_dir / p).resolve()
        bases.append(p)

    subs = _enumerate_subjects(bases)
    print(f"[INFO] find {len(subs)} subjects, for example {subs[:5]}")

    rows = [summarize_subject(sub, bases, TASK, MAX_RUNS) for sub in subs]
    df = pd.DataFrame(rows)

    cols = ["subid"]
    for r in range(1, MAX_RUNS + 1):
        #  summary information
        cols += [
            f"rest{r}_meanFD",
            f"rest{r}_nFrames",
            f"rest{r}_propFD_leq_0p2",
            f"rest{r}_keep",
        ]
    cols += ["kept_runs_count", "meanFD_kept_runs", "check_has2runs"]
    df = df[cols]

    out_path = Path(OUTPUT)
    if not out_path.is_absolute():
        out_path = (script_dir / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format="%.6f")

    # sublist.txt
    eligible = df.loc[df["check_has2runs"] == 1, "subid"].astype(str).tolist()
    sublist_out = out_path.with_name("sublist.txt")
    with sublist_out.open("w", encoding="utf-8") as f:
        for sid in eligible:
            f.write(f"{sid}\n")

    print(f"[OK] Wrote {out_path} with {len(df)} rows.")
    print(f"[OK] Wrote {sublist_out} with {len(eligible)} subjects (>=2 runs kept).")


if __name__ == "__main__":
    main()