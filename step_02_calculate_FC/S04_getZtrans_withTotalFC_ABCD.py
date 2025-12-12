#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import numpy as np
import pandas as pd
from typing import Tuple, Optional, List


CSV_PATH    = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/5th_prediction/prepare_data/cognition/abcd_mergedTable_with_motion2run_deNaN_from_n5959_QCpassed_n4388_label.csv"                     
GM_DIR      = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth/GM_FC"        
WM_DIR      = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth/WM_FC"            
GMWM_DIR    = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth/GM_WM_FC"         
OUT_DIR     = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth_merge"                   
TRI         = "lower"                          
SUBJECT_COL = "subject_id"                   
FILL_VALUE  = 0.0                               



def read_subject_ids(csv_path: str, col: str = "subject_id") -> List[str]:
    df = pd.read_csv(csv_path)
    if col not in df.columns:
        raise ValueError(f"csv not found: {col}")
    return [str(x).strip() for x in df[col].tolist()]

def normalize_sid(s: str) -> str:
    # remove _
    return re.sub(r"[^A-Za-z0-9]", "", s)

def fisher_z(mat: np.ndarray) -> np.ndarray:
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.arctanh(mat)

def replace_nonfinite(arr: np.ndarray, fill: float = 0.0) -> np.ndarray:
    return np.nan_to_num(arr, nan=fill, posinf=fill, neginf=fill)

def vectorize_symmetric(mat: np.ndarray, tri: str = "lower") -> np.ndarray:
    assert mat.ndim == 2 and mat.shape[0] == mat.shape[1]
    n = mat.shape[0]
    idx = np.triu_indices(n, 1) if tri == "upper" else np.tril_indices(n, -1)
    return mat[idx].ravel()

def vectorize_rectangular(mat: np.ndarray) -> np.ndarray:
    return mat.reshape(1, -1).ravel()

def fc_paths(sid_norm: str) -> Tuple[str, str, str]:
    gm = os.path.join(GM_DIR,   f"{sid_norm}_GM_FC.npy")
    wm = os.path.join(WM_DIR,   f"{sid_norm}_WM_FC.npy")
    gw = os.path.join(GMWM_DIR, f"{sid_norm}_GM_WM_FC.npy")
    return gm, wm, gw

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_ids = read_subject_ids(CSV_PATH, SUBJECT_COL)
    print(f"CSV {len(all_ids)} subject_id")

    used_ids: List[str] = []
    gm_vecs, wm_vecs, gw_vecs = [], [], []
    gm_shape: Optional[Tuple[int, int]] = None
    wm_shape: Optional[Tuple[int, int]] = None
    gw_shape: Optional[Tuple[int, int]] = None
    skipped = 0

    for i, sid in enumerate(all_ids, 1):
        sid_norm = normalize_sid(sid)
        gm_p, wm_p, gw_p = fc_paths(sid_norm)

        if not (os.path.exists(gm_p) and os.path.exists(wm_p) and os.path.exists(gw_p)):
            print(f"[skip] {sid} -> not found file:"
                  f"{'GM ' if not os.path.exists(gm_p) else ''}"
                  f"{'WM ' if not os.path.exists(wm_p) else ''}"
                  f"{'GW' if not os.path.exists(gw_p) else ''}".strip())
            skipped += 1
            continue

        try:
            gm = np.load(gm_p); wm = np.load(wm_p); gw = np.load(gw_p)
        except Exception as e:
            print(f"[skip] {sid} -> np.load error: {e}")
            skipped += 1
            continue

        if gm.ndim != 2 or gm.shape[0] != gm.shape[1]:
            print(f"[skip] {sid} -> GM, shape={gm.shape}"); skipped += 1; continue
        if wm.ndim != 2 or wm.shape[0] != wm.shape[1]:
            print(f"[skip] {sid} -> WM, shape={wm.shape}"); skipped += 1; continue
        if gw.ndim != 2:
            print(f"[skip] {sid} -> GW, shape={gw.shape}"); skipped += 1; continue

        gm_z, wm_z, gw_z = fisher_z(gm), fisher_z(wm), fisher_z(gw)

        if gm_shape is None: gm_shape = gm_z.shape
        if wm_shape is None: wm_shape = wm_z.shape
        if gw_shape is None: gw_shape = gw_z.shape

        if gm_z.shape != gm_shape or wm_z.shape != wm_shape or gw_z.shape != gw_shape:
            print(f"[skip] {sid} -> Inconsistent shapes: GM{gm_z.shape}!={gm_shape} or "
                  f"WM{wm_z.shape}!={wm_shape} or GW{gw_z.shape}!={gw_shape}")
            skipped += 1
            continue

        gm_vec = vectorize_symmetric(gm_z, TRI).astype(np.float32)
        wm_vec = vectorize_symmetric(wm_z, TRI).astype(np.float32)
        gw_vec = vectorize_rectangular(gw_z).astype(np.float32)

        gm_vec = replace_nonfinite(gm_vec, FILL_VALUE)
        wm_vec = replace_nonfinite(wm_vec, FILL_VALUE)
        gw_vec = replace_nonfinite(gw_vec, FILL_VALUE)


        gm_vecs.append(gm_vec); wm_vecs.append(wm_vec); gw_vecs.append(gw_vec)
        used_ids.append(sid)

        if i % 50 == 0 or i == len(all_ids):
            print(f"[Progress] Processing {i}/{len(all_ids)}")

    if not used_ids:
        raise RuntimeError("No subjects passed the test.")

    GM = replace_nonfinite(np.vstack(gm_vecs), FILL_VALUE)
    WM = replace_nonfinite(np.vstack(wm_vecs), FILL_VALUE)
    GW = replace_nonfinite(np.vstack(gw_vecs), FILL_VALUE)

    gm_out  = os.path.join(OUT_DIR, "total_FCvector_GG_n4388_cognition.txt")
    wm_out  = os.path.join(OUT_DIR, "total_FCvector_WW_60p_n4388_cognition.txt")
    gw_out  = os.path.join(OUT_DIR, "total_FCvector_GW_60p_n4388_cognition.txt")
    sid_out = os.path.join(OUT_DIR, "subject_ids_used_n4388_cognition_motion2runFD.txt")

    np.savetxt(gm_out, GM, fmt="%.6f", delimiter="\t")
    np.savetxt(wm_out, WM, fmt="%.6f", delimiter="\t")
    np.savetxt(gw_out, GW, fmt="%.6f", delimiter="\t")
    with open(sid_out, "w", encoding="utf-8") as f:
        f.writelines([s + "\n" for s in used_ids])

    print(f"Valid participants:: {len(used_ids)} / {len(all_ids)}; skipping: {skipped}")
    print(f"GM vector shape: {GM.shape} | WM vector shape: {WM.shape} | G-W vector shape: {GW.shape}")
    print(f"output:\n  {gm_out}\n  {wm_out}\n  {gw_out}\n  {sid_out}")

if __name__ == "__main__":
    main()