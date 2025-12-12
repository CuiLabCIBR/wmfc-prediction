#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
from pathlib import Path
from typing import Optional, List, Tuple
import os, re, sys, shutil, subprocess
import numpy as np
import pandas as pd
import nibabel as nib

FMRIPREP_DIRS: List[Path] = [
    Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/fmriprep"),
    # Path("/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets_new/brainproject/fmriprep"),
    # Path("/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/brainproject/fmriprep_combined"),
]
XCPD_DIRS: List[Path] = [
    Path("/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets_new/brainproject/xcpd/step_2nd_24PcsfGlobal"),
    # Path("/ibmgpfs/cuizaixu_lab/xuhaoshu/WM_prediction/datasets/brainproject/xcpd/step_2nd_24PcsfGlobal"),
]

SUBLIST_TSV  = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/sublist_runs_new.tsv")   
OUT_ROOT     = Path("/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/smooth")     

TASK         = "rest"
SPACE        = "MNI152NLin6Asym"
RES          = "2"
FWHM_MM      = 4.0
GM_LABEL     = 1   # 1: GM
WM_LABEL     = 2   # 2: WM


WB_CMD = os.environ.get("WB_CMD", "wb_command")


RUN_RE = re.compile(r"run[-_]?(\d+)", re.IGNORECASE)

def _norm_subid(s: str) -> str:
    s = str(s).strip()
    return s if s.startswith("sub-") else f"sub-{s}"

def _need_wb():
    if shutil.which(WB_CMD) is None:
        print("[ERR] wb_command not found.  WB_CMD=/path/to/wb_command", file=sys.stderr)
        sys.exit(1)

def _find_dseg(subid: str) -> Optional[Path]:
    rel = Path(subid) / "anat" / f"{subid}_run-1_space-{SPACE}_res-{RES}_dseg.nii.gz"
    for root in FMRIPREP_DIRS:
        p = root / rel
        if p.exists():
            return p
    return None

def _find_bold(subid: str, run_token: str) -> Optional[Path]:
    m = RUN_RE.search(run_token)
    run_num = m.group(1) if m else None
    if not run_num:
        return None
    rel = Path(subid) / "func" / f"{subid}_task-{TASK}_run-{run_num}_space-{SPACE}_res-{RES}_desc-denoised_bold.nii.gz"
    for root in XCPD_DIRS:
        p = root / rel
        if p.exists():
            return p
    return None

def _load(path: Path) -> Tuple[nib.Nifti1Image, np.ndarray]:
    img = nib.load(str(path))
    data = img.get_fdata(dtype=np.float32)
    return img, data

def _save_like(ref_img: nib.Nifti1Image, data: np.ndarray, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out = nib.Nifti1Image(data, ref_img.affine, ref_img.header)
    out.set_data_dtype(np.float32 if data.dtype.kind == "f" else data.dtype)
    nib.save(out, str(out_path))

def _run(cmd: list[str]):
    print("[CMD]", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

def process_one(subid: str, run_token: str):
    subid = _norm_subid(subid)
    out_dir = OUT_ROOT / subid
    out_dir.mkdir(parents=True, exist_ok=True)

    dseg = _find_dseg(subid)
    if dseg is None:
        print(f"[WARN] {subid}: dseg not found -> skip")
        return

    bold = _find_bold(subid, run_token)
    if bold is None:
        print(f"[WARN] {subid} {run_token}: denoised BOLD not found -> skip")
        return

    dseg_img, dseg_data = _load(dseg)
    gm_mask = (dseg_data == GM_LABEL)
    wm_mask = (dseg_data == WM_LABEL)


    base_mask = f"{subid}_space-{SPACE}_res-{RES}_dseg"
    gm_mask_path = out_dir / f"{base_mask}_GM.nii.gz"
    wm_mask_path = out_dir / f"{base_mask}_WM.nii.gz"
    if not gm_mask_path.exists():
        _save_like(dseg_img, gm_mask.astype(np.uint8), gm_mask_path)
    if not wm_mask_path.exists():
        _save_like(dseg_img, wm_mask.astype(np.uint8), wm_mask_path)

    bold_img, bold_4d = _load(bold)
    if bold_4d.ndim != 4:
        print(f"[WARN] {subid} {run_token}: BOLD not 4D -> skip")
        return
    if bold_4d.shape[:3] != gm_mask.shape:
        print(f"[WARN] {subid} {run_token}: dseg and BOLD shape mismatch -> skip")
        return

    gm_4d = bold_4d * gm_mask[..., None]    
    wm_4d = bold_4d * wm_mask[..., None]

    stem = bold.name.replace("_desc-denoised_bold.nii.gz", "")
    gm_4d_path = out_dir / f"{stem}_GM.nii.gz"
    wm_4d_path = out_dir / f"{stem}_WM.nii.gz"
    _save_like(bold_img, gm_4d, gm_4d_path)
    _save_like(bold_img, wm_4d, wm_4d_path)
    print("GM and WM bold data getted!")
    print(f"[OK] {subid} {run_token}: Final masked outputs saved -> {out_dir}")

def main():
    _need_wb()
    sublist = Path(os.environ.get("OVERRIDE_SUBLIST", "") or SUBLIST_TSV)
    if not sublist.exists():
        print(f"[ERR] sublist TSV not found: {sublist}", file=sys.stderr); sys.exit(1)

    df = pd.read_csv(sublist, sep=None, engine="python")
    if not {"subid","run"}.issubset(df.columns):
        print("[ERR] sublist_runs.tsv must have columns: subid, run", file=sys.stderr); sys.exit(1)

    for subid, run_token in df[["subid","run"]].itertuples(index=False, name=None):
        subid = str(subid).strip(); run_token = str(run_token).strip()
        if subid and run_token:
            process_one(subid, run_token)

if __name__ == "__main__":
    main()
