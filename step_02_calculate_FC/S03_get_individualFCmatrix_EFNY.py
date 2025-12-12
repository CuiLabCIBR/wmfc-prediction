import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib

DATA_ROOT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/smooth_CSF"
OUT_ROOT  = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/FC_nosmooth_CSF"

SUBJECTS = []
SUBJECTS_FILE = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/results/sublist_final_532.txt"

RUNS_CAND = ["1", "2", "3", "4"]

NAMING = {"task": "rest", "space": "MNI152NLin6Asym"}

ATLAS_GM_PATH = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/3rd_cal_FC/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm_resample.nii.gz"
ATLAS_WM_PATH = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/EFNY/code/3rd_cal_FC/rICBM_DTI_81_WMPM_60p_FMRIB58_resample.nii.gz"

# GM_FILE_TEMPLATE = "{sub}/{sub}_task-{task}_run-{run}_space-{space}_res-2_GM_smoothMaksed.nii.gz"
# WM_FILE_TEMPLATE = "{sub}/{sub}_task-{task}_run-{run}_space-{space}_res-2_WM_smoothMaksed.nii.gz"
GM_FILE_TEMPLATE = "{sub}/{sub}_task-{task}_run-{run}_space-{space}_res-2_GM.nii.gz"
WM_FILE_TEMPLATE = "{sub}/{sub}_task-{task}_run-{run}_space-{space}_res-2_WM.nii.gz"

SAVE_TSV = True
SAVE_NPY = True

RUN_USAGE_TSV = os.path.join(OUT_ROOT, "run_usage_summary.tsv")       
RUN_USAGE_PARTS_DIR = os.path.join(OUT_ROOT, "run_usage_parts")          

def ensure_dir(d):
    os.makedirs(d, exist_ok=True); return d

def list_subjects():
    one = os.environ.get("SUBJECT", "").strip()
    if one: return [one]
    if SUBJECTS_FILE and os.path.exists(SUBJECTS_FILE):
        with open(SUBJECTS_FILE, "r", encoding="utf-8") as f:
            subs = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
        if subs: return subs
    if SUBJECTS: return SUBJECTS
    subs = [d for d in os.listdir(DATA_ROOT) if os.path.isdir(os.path.join(DATA_ROOT, d))]
    subs.sort(); return subs

def path_if_exists(template, sub, run):
    rel = template.format(sub=sub, run=run, **NAMING)
    p = os.path.join(DATA_ROOT, rel)
    return p if os.path.exists(p) else None

def collect_available_runs(sub):
    gm_list, wm_list, used_runs, skipped_info = [], [], [], []
    for r in RUNS_CAND:
        g = path_if_exists(GM_FILE_TEMPLATE, sub, r)
        w = path_if_exists(WM_FILE_TEMPLATE, sub, r)
        if g and w:
            gm_list.append(g); wm_list.append(w); used_runs.append(r)
        else:
            if g and not w:
                skipped_info.append(f"{r}:only_GM")
                print(f"[WARN] {sub} run-{r}: Only GM found, this run is skipped.")
            elif w and not g:
                skipped_info.append(f"{r}:only_WM")
                print(f"[WARN] {sub} run-{r}: Only WM found, this run is skipped.")
            else:
                skipped_info.append(f"{r}:missing")
    return gm_list, wm_list, used_runs, skipped_info

def load_concat_4d(paths):
    imgs = [nib.load(p) for p in paths]
    shapes = [img.shape[:3] for img in imgs]
    if len(set(shapes)) != 1:
        raise ValueError(f"The space dimensions vary between different runs:{shapes}")
    affines = [img.affine for img in imgs]
    for a in affines[1:]:
        if not np.allclose(affines[0], a, atol=1e-3):
            raise ValueError("The affine matrices for different runs are inconsistent.")
    arrays = [img.get_fdata(dtype=np.float32) for img in imgs]
    data4d = np.concatenate(arrays, axis=3)
    return data4d, affines[0]

def check_same_space(data_shape_xyz, data_affine, atlas_img, tag):
    if data_shape_xyz != atlas_img.shape[:3]:
        raise ValueError(f"[{tag}] atlas voxel {atlas_img.shape[:3]} != Functional data {data_shape_xyz}")
    if not np.allclose(data_affine, atlas_img.affine, atol=1e-3):
        raise ValueError(f"[{tag}] atlas and functional data affine mapping are inconsistent")

def atlas_roi_indices(atlas_img):
    atl = atlas_img.get_fdata()
    labels = np.unique(atl)
    labels = labels[np.isfinite(labels)].astype(int)
    labels = labels[labels != 0]
    labels.sort()
    flat = atl.reshape(-1)
    idx_list = [np.where(flat == lab)[0] for lab in labels]
    return labels, idx_list

def roi_mean_timeseries(data4d, idx_list):
    T = data4d.shape[3]
    V = np.prod(data4d.shape[:3])
    data2d = data4d.reshape(V, T).T  # (T, V)
    N = len(idx_list)
    ts = np.empty((T, N), dtype=np.float32)
    for i, idx in enumerate(idx_list):
        ts[:, i] = np.nan if idx.size == 0 else data2d[:, idx].mean(axis=1)
    return ts

def corrcoef_blocks(gm_ts, wm_ts):
    Ngm, Nwm = gm_ts.shape[1], wm_ts.shape[1]
    all_ts = np.hstack([gm_ts, wm_ts]).astype(np.float64, copy=False)
    C = np.corrcoef(all_ts, rowvar=False)
    gm_fc   = np.where(np.isfinite(C[:Ngm, :Ngm]),   np.clip(C[:Ngm, :Ngm],   -1.0, 1.0), np.nan)
    wm_fc   = np.where(np.isfinite(C[Ngm:, Ngm:]),   np.clip(C[Ngm:, Ngm:],   -1.0, 1.0), np.nan)
    gmwm_fc = np.where(np.isfinite(C[:Ngm, Ngm:]),   np.clip(C[:Ngm, Ngm:],   -1.0, 1.0), np.nan)
    return gm_fc, wm_fc, gmwm_fc

def write_run_usage_rows(rows, out_path):
    df = pd.DataFrame(rows, columns=["subid", "n_runs", "runs", "skipped"])
    df.to_csv(out_path, sep="\t", index=False)

def write_run_usage_part(row, parts_dir):
    ensure_dir(parts_dir)
    out_path = os.path.join(parts_dir, f"{row['subid']}.tsv")
    pd.DataFrame([row], columns=["subid","n_runs","runs","skipped"]).to_csv(out_path, sep="\t", index=False)
    print(f"[STATS] Written: {out_path}")

def merge_run_usage_parts(parts_dir, out_path):
    files = sorted(glob.glob(os.path.join(parts_dir, "*.tsv")))
    if not files:
        print(f"[STATS] not found:{parts_dir}")
        return
    dfs = [pd.read_csv(f, sep="\t") for f in files]
    df = pd.concat(dfs, ignore_index=True)

    df = df.drop_duplicates(subset=["subid"], keep="last").sort_values(["subid"]).reset_index(drop=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[STATS] has been merged {len(files)} --> {out_path}")

def main():
    out_gm    = ensure_dir(os.path.join(OUT_ROOT, "GM_FC"))
    out_wm    = ensure_dir(os.path.join(OUT_ROOT, "WM_FC"))
    out_gmwm  = ensure_dir(os.path.join(OUT_ROOT, "GM_WM_FC"))
    out_ts_gm = ensure_dir(os.path.join(OUT_ROOT, "GM_TS"))
    out_ts_wm = ensure_dir(os.path.join(OUT_ROOT, "WM_TS"))
    ensure_dir(OUT_ROOT)
    ensure_dir(RUN_USAGE_PARTS_DIR)

    if os.environ.get("MERGE_RUN_USAGE", "") == "1":
        merge_run_usage_parts(RUN_USAGE_PARTS_DIR, RUN_USAGE_TSV)
        return

    atlas_gm = nib.load(ATLAS_GM_PATH)
    atlas_wm = nib.load(ATLAS_WM_PATH)
    gm_labels, gm_idx = atlas_roi_indices(atlas_gm)
    wm_labels, wm_idx = atlas_roi_indices(atlas_wm)
    gm_names = [f"GM_{int(x)}" for x in gm_labels]
    wm_names = [f"WM_{int(x)}" for x in wm_labels]

    subs = list_subjects()
    print(f"processing {len(subs)} subjects: {subs}")

    run_usage_rows = []
    for sub in subs:
        print(f"\n==== {sub} ====")
        gm_paths, wm_paths, used_runs, skipped_info = collect_available_runs(sub)

        # record the row
        row = {
            "subid": sub,
            "n_runs": len(used_runs),
            "runs": ",".join(used_runs) if used_runs else "",
            "skipped": ";".join(skipped_info) if skipped_info else ""
        }
        run_usage_rows.append(row)

        if not gm_paths:
            print(f"[WARN] {sub}: No run available (GM/WM pairs exist) -> Skip")
        
            if os.environ.get("SUBJECT","").strip():
                write_run_usage_part(row, RUN_USAGE_PARTS_DIR)
            continue

        # ====== calculation ======
        gm_4d, gm_aff = load_concat_4d(gm_paths)
        wm_4d, wm_aff = load_concat_4d(wm_paths)
        check_same_space(gm_4d.shape[:3], gm_aff, atlas_gm, "GM")
        check_same_space(wm_4d.shape[:3], wm_aff, atlas_wm, "WM")
        gm_ts = roi_mean_timeseries(gm_4d, gm_idx)
        wm_ts = roi_mean_timeseries(wm_4d, wm_idx)
        if gm_ts.shape[0] != wm_ts.shape[0]:
            raise ValueError(f"T Inconsistency: GM={gm_ts.shape[0]} vs WM={wm_ts.shape[0]}")
        gm_fc, wm_fc, gmwm_fc = corrcoef_blocks(gm_ts, wm_ts)

        # ====== save ======
        np.save(os.path.join(out_ts_gm, f"{sub}_GM_timeseries.npy"), gm_ts) if SAVE_NPY else None
        np.save(os.path.join(out_ts_wm, f"{sub}_WM_timeseries.npy"), wm_ts) if SAVE_NPY else None
        if SAVE_TSV:
            pd.DataFrame(gm_ts, index=[f"t{i}" for i in range(gm_ts.shape[0])], columns=gm_names).to_csv(
                os.path.join(out_ts_gm, f"{sub}_GM_timeseries.tsv"), sep="\t", float_format="%.6f", na_rep="NaN")
            pd.DataFrame(wm_ts, index=[f"t{i}" for i in range(wm_ts.shape[0])], columns=wm_names).to_csv(
                os.path.join(out_ts_wm, f"{sub}_WM_timeseries.tsv"), sep="\t", float_format="%.6f", na_rep="NaN")

        def save_matrix(mat, base, rows, cols):
            if SAVE_NPY: np.save(base + ".npy", mat)
            if SAVE_TSV: pd.DataFrame(mat, index=rows, columns=cols).to_csv(base + ".tsv", sep="\t", float_format="%.6f", na_rep="NaN")

        save_matrix(gm_fc,   os.path.join(out_gm,   f"{sub}_GM_FC"),    gm_names, gm_names)
        save_matrix(wm_fc,   os.path.join(out_wm,   f"{sub}_WM_FC"),    wm_names, wm_names)
        save_matrix(gmwm_fc, os.path.join(out_gmwm, f"{sub}_GM_WM_FC"), gm_names, wm_names)
        print("complete")

        
        if os.environ.get("SUBJECT","").strip():
            write_run_usage_part(row, RUN_USAGE_PARTS_DIR)

    # summary FC runs
    if not os.environ.get("SUBJECT","").strip():
        write_run_usage_rows(run_usage_rows, RUN_USAGE_TSV)
        print(f"\n[STATS]  run have been saved: {RUN_USAGE_TSV}")
    else:
        print(f"\n[STATS]  {RUN_USAGE_PARTS_DIR}")

if __name__ == "__main__":
    main()
