#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
import re
import sys

# ========= path =========
CSV_PATH     = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/head_motion/fd_summary_addframe.csv")
PAIRS_OUT    = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/sublist_runs_new.tsv")
SUBLIST_PATH = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/batchcode/resubmit.txt")
# =================================


def detect_run_indices(columns):
    """ restX_keep """
    run_ids = []
    pat = re.compile(r"^rest(\d+)_keep$")
    for col in columns:
        m = pat.match(col)
        if m:
            run_ids.append(int(m.group(1)))
    run_ids.sort()
    return run_ids


def read_sublist(path: Path):
    """read sublist.txt  ID"""
    subs = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        subs.add(s)
    return subs


def _to_bool_keep(val) -> bool:

    if pd.isna(val):
        return False
    try:
        return float(val) == 1.0
    except Exception:
        return str(val).strip() in {"1", "True", "true", "TRUE"}


def _frames_equal_179(val) -> bool:
    """Determine if nFrames is 179."""
    if pd.isna(val):
        return False
    try:
        return int(float(val)) == 179
    except Exception:
        return False


def main():
    if not CSV_PATH.exists():
        print(f"[error] could not find file:{CSV_PATH}", file=sys.stderr)
        sys.exit(1)

    # read CSV
    try:
        df = pd.read_csv(CSV_PATH, sep=None, engine="python")
    except Exception as e:
        print(f"Failed to read CSV:{e}", file=sys.stderr)
        sys.exit(1)

    if "subid" not in df.columns:
        print("[error] Required columns are missing:'subid'", file=sys.stderr)
        sys.exit(1)


    run_ids = detect_run_indices(df.columns)
    if not run_ids:
        print("[error] No column named 'restX_keep' was found, such as 'rest1_keep'", file=sys.stderr)
        sys.exit(1)

    # according to SUBLIST_PATH or check_has2runs select subjects
    if SUBLIST_PATH is not None:
        if not Path(SUBLIST_PATH).exists():
            print(f"[error] SUBLIST_PATH does not exist: {SUBLIST_PATH}", file=sys.stderr)
            sys.exit(1)
        keep_set = read_sublist(Path(SUBLIST_PATH))
        df_sel = df[df["subid"].astype(str).isin(keep_set)].copy()
        filter_info = f"SUBLIST_PATH Selected participants: {df_sel['subid'].nunique()}"
    else:
        if "check_has2runs" not in df.columns:
            print("[error] SUBLIST_PATH is not provided and the 'check_has2runs' column is missing.", file=sys.stderr)
            sys.exit(1)
        df_sel = df[df["check_has2runs"] == 1].copy()
        filter_info = f"Selected participants by using check_has2runs==1:{df_sel['subid'].nunique()}"

    # ===== Further filter by "valid run" (keep==1 and nFrames==179), requiring >= 2 valid runs. =====
    pairs = []
    n_valid_subjects = 0

    for _, row in df_sel.iterrows():
        sub = str(row["subid"])
        valid_runs_for_sub = []

        for i in run_ids:
            col_keep    = f"rest{i}_keep"
            col_nframes = f"rest{i}_nFrames"

            if col_keep not in df_sel.columns or col_nframes not in df_sel.columns:
                continue

            keep_flag   = _to_bool_keep(row[col_keep])
            frames_flag = _frames_equal_179(row[col_nframes])

            # Effective run: keep==1 and nFrames==179
            if keep_flag and frames_flag:
                valid_runs_for_sub.append(i)

        # Only participants with a valid run >= 2 were retained.
        if len(valid_runs_for_sub) >= 2:
            n_valid_subjects += 1
            for i in valid_runs_for_sub:
                pairs.append({"subid": sub, "run": f"rest_run-{i}"})

    if pairs:
        out_df = pd.DataFrame(pairs).sort_values(["subid", "run"])
    else:
        out_df = pd.DataFrame(columns=["subid", "run"])

    # save (subid, run) mapping table
    out_df.to_csv(PAIRS_OUT, sep="\t", index=False)

    # ===== Save the final subids that meet the criteria to a txt file. =====
    if not out_df.empty:

        unique_subids = sorted(out_df["subid"].astype(str).unique())

        sublist_out = PAIRS_OUT.with_name("sublist_final_533.txt")

        with sublist_out.open("w", encoding="utf-8") as f:
            for sid in unique_subids:
                f.write(f"{sid}\n")

        print(f"[Output] Final results that meet the conditions:{sublist_out.resolve()}")
    else:
        print("[WARNING] No subjects met the criteria; no sublist txt file was generated.")

    print(f"[Completed] {filter_info}")
    print(f"[Completed] Number of subjects satisfying "valid run >= 2 and each valid run has nFrames = 179 and keep = 1":{n_valid_subjects}")
    print(f"[Completed] Number of qualified (subid, run) records:{len(pairs)}")
    print(f"[Output] Qualified run mapping:{PAIRS_OUT.resolve()}")


if __name__ == "__main__":
    main()