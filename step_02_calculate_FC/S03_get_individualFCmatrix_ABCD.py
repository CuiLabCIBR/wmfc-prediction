import os
import numpy as np
import pandas as pd
import nibabel as nib

DATA_ROOT = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/smooth"          
OUT_ROOT  = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/FC_nosmooth"  


SUBJECTS = []              # subid/[] []:list all sub
SUBJECTS_FILE = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/check_data/has_smooth_data_new.txt"       #file "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/results/check_data/has_smooth_data_new.txt"/None

RUNS = ["1", "2"]

# naming
NAMING = {
    "ses":   "baselineYear1Arm1",
    "task":  "rest",
    "space": "MNI152NLin6Asym",
    "desc":  "denoised",
}

# atlas
ATLAS_GM_PATH = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/4th_calFC/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm_resample.nii.gz"
ATLAS_WM_PATH = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/4th_calFC/rICBM_DTI_81_WMPM_60p_FMRIB58_resample.nii.gz"

# GM_FILE_TEMPLATE = "sub-{sub}/sub-{sub}_ses-{ses}_task-{task}_run-{run}_space-{space}_desc-{desc}_GM_smoothMasked.nii.gz"
# WM_FILE_TEMPLATE = "sub-{sub}/sub-{sub}_ses-{ses}_task-{task}_run-{run}_space-{space}_desc-{desc}_WM_smoothMasked.nii.gz"
GM_FILE_TEMPLATE = "sub-{sub}/sub-{sub}_ses-{ses}_task-{task}_run-{run}_space-{space}_desc-{desc}_GM.nii.gz"
WM_FILE_TEMPLATE = "sub-{sub}/sub-{sub}_ses-{ses}_task-{task}_run-{run}_space-{space}_desc-{desc}_WM.nii.gz"

# out
SAVE_TSV = True
SAVE_NPY = True

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)
    return d


def list_subjects():

    env_one = os.environ.get("SUBJECT", "").strip()
    if env_one:
        return [env_one]

    if SUBJECTS_FILE and os.path.exists(SUBJECTS_FILE):
        with open(SUBJECTS_FILE, "r", encoding="utf-8") as f:
            subs = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
        if subs:
            return subs

    if SUBJECTS:
        return SUBJECTS

    subs = [d for d in os.listdir(DATA_ROOT) if os.path.isdir(os.path.join(DATA_ROOT, d))]
    subs.sort()
    return subs

def build_path(template, sub, run):
    rel = template.format(sub=sub, run=run, **NAMING)
    path = os.path.join(DATA_ROOT, rel)
    if not os.path.exists(path):
        raise FileNotFoundError(f"not found：{path}")
    return path

def load_concat_4d(paths):
    imgs = [nib.load(p) for p in paths]
    shapes = [img.shape[:3] for img in imgs]
    if len(set(shapes)) != 1:
        raise ValueError(f"run shape different：{shapes}")
    affines = [img.affine for img in imgs]
    for a in affines[1:]:
        if not np.allclose(affines[0], a, atol=1e-3):
            raise ValueError("run affine different")
    arrays = [img.get_fdata(dtype=np.float32) for img in imgs]  
    data4d = np.concatenate(arrays, axis=3)               #concatenate   
    return data4d, affines[0]

def check_same_space(data_shape_xyz, data_affine, atlas_img, tag):
    if data_shape_xyz != atlas_img.shape[:3]:
        raise ValueError(f"[{tag}] atlas voxel size {atlas_img.shape[:3]} != func {data_shape_xyz}")
    if not np.allclose(data_affine, atlas_img.affine, atol=1e-3):
        raise ValueError(f"[{tag}] atlas different from func affine")

def atlas_roi_indices(atlas_img):
    # (labels_sorted, idx_list)
    atl = atlas_img.get_fdata()
    labels = np.unique(atl)
    labels = labels[np.isfinite(labels)].astype(int)
    labels = labels[labels != 0]
    labels.sort()
    flat = atl.reshape(-1)
    idx_list = [np.where(flat == lab)[0] for lab in labels]
    return labels, idx_list

def roi_mean_timeseries(data4d, idx_list):

    # data4d: (X,Y,Z,T)     idx_list
    # T x Nroi：timeseries

    T = data4d.shape[3]
    V = np.prod(data4d.shape[:3])
    data2d = data4d.reshape(V, T).T  # (T, V)
    N = len(idx_list)
    ts = np.empty((T, N), dtype=np.float32)
    for i, idx in enumerate(idx_list):
        ts[:, i] = np.nan if idx.size == 0 else data2d[:, idx].mean(axis=1)
    return ts

def corrcoef_blocks(gm_ts, wm_ts):
  
    #Pearson ：gm_ts: (T x Ngm), wm_ts: (T x Nwm)
    #gm_fc (Ngm x Ngm), wm_fc (Nwm x Nwm), gmwm_fc (Ngm x Nwm)

    Ngm, Nwm = gm_ts.shape[1], wm_ts.shape[1]
    all_ts = np.hstack([gm_ts, wm_ts]).astype(np.float64, copy=False)  # (T, Ngm+Nwm)
    C = np.corrcoef(all_ts, rowvar=False)                               # (Ngm+Nwm)×(Ngm+Nwm)

    gm_fc   = C[:Ngm, :Ngm]
    wm_fc   = C[Ngm:, Ngm:]
    gmwm_fc = C[:Ngm, Ngm:]

    gm_fc   = np.where(np.isfinite(gm_fc),   np.clip(gm_fc,   -1.0, 1.0), np.nan)
    wm_fc   = np.where(np.isfinite(wm_fc),   np.clip(wm_fc,   -1.0, 1.0), np.nan)
    gmwm_fc = np.where(np.isfinite(gmwm_fc), np.clip(gmwm_fc, -1.0, 1.0), np.nan)
    return gm_fc, wm_fc, gmwm_fc



def save_matrix(mat, base, row_labels=None, col_labels=None):
    if SAVE_NPY:
        np.save(base + ".npy", mat)
    if SAVE_TSV:
        idx = row_labels if row_labels is not None else [f"r{i}" for i in range(mat.shape[0])]
        cols = col_labels if col_labels is not None else [f"c{j}" for j in range(mat.shape[1])]
        df = pd.DataFrame(mat, index=idx, columns=cols)
        df.to_csv(base + ".tsv", sep="\t", float_format="%.6f", na_rep="NaN")

def save_timeseries(ts, base, col_labels):
    if SAVE_NPY:
        np.save(base + ".npy", ts)
    if SAVE_TSV:
        T = ts.shape[0]
        idx = [f"t{i}" for i in range(T)]
        df = pd.DataFrame(ts, index=idx, columns=col_labels)
        df.to_csv(base + ".tsv", sep="\t", float_format="%.6f", na_rep="NaN")
        

def main():
    out_gm   = ensure_dir(os.path.join(OUT_ROOT, "GM_FC"))
    out_wm   = ensure_dir(os.path.join(OUT_ROOT, "WM_FC"))
    out_gmwm = ensure_dir(os.path.join(OUT_ROOT, "GM_WM_FC"))
    out_ts_gm = ensure_dir(os.path.join(OUT_ROOT, "GM_TS"))  
    out_ts_wm = ensure_dir(os.path.join(OUT_ROOT, "WM_TS"))  

    # atlas
    atlas_gm = nib.load(ATLAS_GM_PATH)
    atlas_wm = nib.load(ATLAS_WM_PATH)
    gm_labels, gm_idx = atlas_roi_indices(atlas_gm)
    wm_labels, wm_idx = atlas_roi_indices(atlas_wm)
    gm_names = [f"GM_{int(x)}" for x in gm_labels]
    wm_names = [f"WM_{int(x)}" for x in wm_labels]

    subs = list_subjects()
    print(f"proprecessing {len(subs)} subjects:{subs}")

    for sub in subs:
        print(f"\n==== {sub} ====")

        gm_paths = [build_path(GM_FILE_TEMPLATE, sub, r) for r in RUNS]
        wm_paths = [build_path(WM_FILE_TEMPLATE, sub, r) for r in RUNS]
        gm_4d, gm_aff = load_concat_4d(gm_paths)
        wm_4d, wm_aff = load_concat_4d(wm_paths)

        # check space
        check_same_space(gm_4d.shape[:3], gm_aff, atlas_gm, "GM")
        check_same_space(wm_4d.shape[:3], wm_aff, atlas_wm, "WM")

        # （T×Nroi）
        gm_ts = roi_mean_timeseries(gm_4d, gm_idx)  # (T, Ngm)
        wm_ts = roi_mean_timeseries(wm_4d, wm_idx)  # (T, Nwm)

        # timeseries
        if gm_ts.shape[0] != wm_ts.shape[0]:
            raise ValueError(f"T GM={gm_ts.shape[0]} vs WM={wm_ts.shape[0]}")

        # corr
        gm_fc, wm_fc, gmwm_fc = corrcoef_blocks(gm_ts, wm_ts)

        # save
        save_timeseries(gm_ts, os.path.join(out_ts_gm, f"{sub}_GM_timeseries"), gm_names)
        save_timeseries(wm_ts, os.path.join(out_ts_wm, f"{sub}_WM_timeseries"), wm_names)
        save_matrix(gm_fc,   os.path.join(out_gm,   f"{sub}_GM_FC"),   gm_names, gm_names)
        save_matrix(wm_fc,   os.path.join(out_wm,   f"{sub}_WM_FC"),   wm_names, wm_names)
        save_matrix(gmwm_fc, os.path.join(out_gmwm, f"{sub}_GM_WM_FC"), gm_names, wm_names)
        print("complete")

    print("\nAll done!")

if __name__ == "__main__":
    main()