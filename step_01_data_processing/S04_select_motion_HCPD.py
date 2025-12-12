#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple

SUBLIST = "/GPFS/cuizaixu_lab_permanent/congjing/hcp_hcpd/batchcode_stable/sub_hcpd.txt"

# fMRIPrep dir
BASES = [
    "/GPFS/cuizaixu_lab_permanent/congjing/hcp_hcpd/results/hcpd_stable_xcpd/xcpd_rest_cifti",
]


OUTPUT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/results/head_motion/fd_summary.csv"


TASK_PREFIX = "rest"  
MAX_RUNS = 4


MOTION_GLOB = "*_motion.tsv"

# 解析：task-REST(\d+) 与 acq-(AP|PA)
# 例：sub-XX_task-REST1_acq-AP_..._motion.tsv
TASK_NUM_PATTERN = re.compile(r"task-{}\s*(\d+)?".format(re.escape(TASK_PREFIX.lower())),
                              re.IGNORECASE)
ACQ_PATTERN = re.compile(r"acq[-_]?([A-Za-z]+)", re.IGNORECASE)

# 将 (REST编号, ACQ方向) 映射到固定槽位 1..4
# 1: REST1_AP, 2: REST1_PA, 3: REST2_AP, 4: REST2_PA
SLOT_MAP = {
    (1, "AP"): 1,
    (1, "PA"): 2,
    (2, "AP"): 3,
    (2, "PA"): 4,
}

def _normalize_subid(raw: str) -> Tuple[str, str]:
    """Return (as_is, with_sub_prefix) forms for searching."""
    s = raw.strip()
    if not s:
        return ("", "")
    if s.startswith("sub-"):
        return (s, s)  # both the same
    return (s, f"sub-{s}")

def _find_confounds_for_sub(
    subid: str,
    bases: List[Path],
    task_prefix: str,
    max_runs: int,
) -> Dict[int, Path]:
    """
    按命名规则为该被试找到最多 max_runs 个 *_motion.tsv，
    并映射到固定槽位（1..4）：
      1→REST1_AP, 2→REST1_PA, 3→REST2_AP, 4→REST2_PA
    """
    as_is, with_sub = _normalize_subid(subid)
    subj_variants = [p for p in [as_is, with_sub] if p]
    found: Dict[int, Path] = {}
    no_run_bucket: List[Path] = []

    for base in bases:
        for subj in subj_variants:
            root = (base / subj)
            if not root.exists():
                continue

            for path in root.rglob(MOTION_GLOB):
                # 只接受 func 目录下的文件
                if "func" not in [p.name for p in path.parents]:
                    continue

                lower_name = path.name.lower()

                # 文件名中必须含 subID（避免串号）
                if as_is and (as_is.lower() not in lower_name) and (with_sub.lower() not in lower_name):
                    continue

                # 匹配 task-REST(\d+)
                m_task = TASK_NUM_PATTERN.search(lower_name)
                if not m_task:
                    continue
                task_num_str = m_task.group(1)
                task_num = int(task_num_str) if task_num_str and task_num_str.isdigit() else 1

                # 匹配 acq-(AP|PA)
                m_acq = ACQ_PATTERN.search(lower_name)
                if not m_acq:
                    # 没有 acq 信息，放入备用桶
                    no_run_bucket.append(path)
                    continue
                acq_dir = m_acq.group(1).upper()
                if acq_dir not in ("AP", "PA"):
                    # 只接受 AP/PA
                    continue

                slot = SLOT_MAP.get((task_num, acq_dir))
                if slot is None:
                    # 未在映射表中（比如 REST3），丢到备用桶
                    no_run_bucket.append(path)
                    continue

                if 1 <= slot <= max_runs and slot not in found:
                    found[slot] = path

            if len(found) >= max_runs:
                break
        if len(found) >= max_runs:
            break

    # 如果仍有空位，用无明确slot的信息补全
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
    """
    从 *_motion.tsv 读取 framewise_displacement 列，返回：
    (mean_fd, prop_leq_0p2, n_frames_used)
    """
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
    task_prefix: str,
    max_runs: int,
):
    runs = _find_confounds_for_sub(subid, bases, task_prefix, max_runs)
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
                # 你原来的保留标准
                if (m <= 0.5) and (p >= 0.40):
                    keep = 1
        result[f"rest{r}_meanFD"] = mean_fd
        result[f"rest{r}_propFD_leq_0p2"] = prop_leq
        result[f"rest{r}_keep"] = keep
        keeps.append(keep)
    result["check_has2runs"] = 1 if sum(keeps) >= 2 else 0
    return result

def main():
    script_dir = Path(__file__).resolve().parent

    sublist_path = Path(SUBLIST)
    if not sublist_path.is_absolute():
        sublist_path = (script_dir / sublist_path).resolve()

    bases = []
    for b in BASES:
        p = Path(b)
        if not p.is_absolute():
            p = (script_dir / p).resolve()
        bases.append(p)

    subs = []
    with sublist_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            subs.append(line)

    rows = [summarize_subject(sub, bases, TASK_PREFIX, MAX_RUNS) for sub in subs]

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

    out_path = Path(OUTPUT)
    if not out_path.is_absolute():
        out_path = (script_dir / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format="%.6f")

    print(f"[OK] Wrote {out_path} with {len(df)} rows.")
    print("Slot mapping: rest1=REST1_AP, rest2=REST1_PA, rest3=REST2_AP, rest4=REST2_PA")

if __name__ == "__main__":
    main()
