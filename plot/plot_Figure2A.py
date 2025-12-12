#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt

# =============================
# datapath
# =============================
pvalue_folder = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/post_prediction/get_significance/Age"
file_folder = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/post_prediction/age_boxplot"
file_path   = os.path.join(file_folder, "age_dataCell_totalSets_forBoxplot.mat")
pvalue_path = os.path.join(pvalue_folder, "Pvalue_byPermutation.mat") 

figure_folder = "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/Figure/Age_raincloud"
os.makedirs(figure_folder, exist_ok=True)
filename1 = "age_4datasets_raincloud_addsig.tiff"
filename2 = "age_4datasets_raincloud_addsig.pdf"


def cell_to_str(x):
    arr = np.array(x).squeeze()
    try:
        if arr.dtype.kind in {"U", "S"} and arr.size == 1:
            return str(arr.item())
        return "".join(arr.astype(str).ravel().tolist())
    except Exception:
        return str(arr)

def cell_to_num1d(x):
    arr = np.array(x).squeeze()
    if arr.size == 0:
        return np.array([], dtype=float)
    try:
        return arr.astype(float).ravel()
    except Exception:
        return np.array([float(el) for el in arr.ravel()])

def cell_to_scalar(x):
    try:
        arr = np.array(x).squeeze()
        return float(arr)
    except Exception:
        return np.nan

# =============================
# read dataCell_totalSets
# =============================
mat = scipy.io.loadmat(file_path)

def find_cell_5x6(md):
    for k in ("dataCell_totalSets", "dataCell.totalSets"):
        if k in md and isinstance(md[k], np.ndarray) and md[k].dtype == object and md[k].shape == (5, 6):
            return md[k]
    for v in md.values():
        if isinstance(v, np.ndarray) and v.dtype == object and v.ndim == 2 and v.shape == (5, 6):
            return v
    return None

data = find_cell_5x6(mat)
if data is None:
    raise RuntimeError("No 5×6 cell found in .mat file (dataCell_totalSets / dataCell.totalSets).")

# =============================
# header 
# =============================
headers = [cell_to_str(data[0, j]) for j in range(data.shape[1])]
cond_map = {
    "R_gg_total":        "G-G",
    "R_gw_total":        "G-W",
    "partialR_gw_total": "G-W/G-G",
    "R_ww_total":        "W-W",
    "partialR_ww_total": "W-W/G-G",
}
cond_order    = ["G-G", "G-W", "G-W/G-G", "W-W", "W-W/G-G"]
dataset_order = ["EFNY", "CCNP", "HCPD", "PNC"]  

# {dataset: {cond: values}}
by_ds = {ds: {c: None for c in cond_order} for ds in dataset_order}
for i in range(1, data.shape[0]):  # row 2..5
    ds = cell_to_str(data[i, 0])
    if ds not in dataset_order:
        continue
    for j in range(1, data.shape[1]):  # column 2..6
        raw = headers[j]
        if raw not in cond_map:
            continue
        cond = cond_map[raw]
        by_ds[ds][cond] = cell_to_num1d(data[i, j])

# DataFrame 表
dfs = []
for ds in dataset_order:
    cols = {}
    for c in cond_order:
        v = by_ds[ds][c]
        if v is None:
            v = np.array([], dtype=float)
        cols[c] = v
    maxlen = max(len(cols[c]) for c in cond_order)
    for c in cond_order:
        if len(cols[c]) < maxlen:
            cols[c] = np.concatenate([cols[c], np.full(maxlen - len(cols[c]), np.nan)])
    dfs.append(pd.DataFrame(cols))

name_list = dataset_order[:]  # x-axis label

# =============================
# read p / q, and generate * map
# =============================
pv = scipy.io.loadmat(pvalue_path)

P_cell = pv.get("P_total", pv.get("p_total", None))
Q_cell = pv.get("Q_total", pv.get("q_total", None))
if P_cell is None or Q_cell is None:
    raise RuntimeError("The variables P_total / Q_total were not found in Pvalue_byPermutation.mat.")

p_headers = [cell_to_str(P_cell[0, j]) for j in range(P_cell.shape[1])]
q_headers = [cell_to_str(Q_cell[0, j]) for j in range(Q_cell.shape[1])]

def idx_in_header(headers, target):
    for i, h in enumerate(headers):
        if h == target:
            return i
    raise ValueError(f"Column names not found in header:{target}")

p_col_idx = {
    "G-G":        idx_in_header(p_headers, "P_GG"),
    "G-W":        idx_in_header(p_headers, "P_GW"),
    "W-W":        idx_in_header(p_headers, "P_WW"),
    "G-W/G-G":    idx_in_header(p_headers, "partialP_GW"),
    "W-W/G-G":    idx_in_header(p_headers, "partialP_WW"),
}
q_col_idx = {
    "G-G":        idx_in_header(q_headers, "q_GG"),
    "G-W":        idx_in_header(q_headers, "q_GW"),
    "W-W":        idx_in_header(q_headers, "q_WW"),
    "G-W/G-G":    idx_in_header(q_headers, "q_partial_GW"),
    "W-W/G-G":    idx_in_header(q_headers, "q_partial_WW"),
}

# row index
row_map = {}
for i in range(1, P_cell.shape[0]):
    ds_name = cell_to_str(P_cell[i, 0])
    row_map[ds_name] = i

# stars_map[dataset][cond] = '','*','**','***'
stars_map = {ds: {c: "" for c in cond_order} for ds in dataset_order}
for ds in dataset_order:
    if ds not in row_map:
        continue
    r = row_map[ds]
    for cond in cond_order:
        p_val = cell_to_scalar(P_cell[r, p_col_idx[cond]])
        q_val = cell_to_scalar(Q_cell[r, q_col_idx[cond]])

        stars = ""
        if not np.isnan(q_val) and q_val <= 0.005:
            stars = "**"
        elif not np.isnan(q_val) and q_val <= 0.05:
            stars = "*"
        stars_map[ds][cond] = stars


# First, determine the y-axis range (based on all values).
all_values_list = []
for df in dfs:
    for c in df.columns:
        vals = pd.to_numeric(df[c], errors='coerce').values
        vals = vals[~np.isnan(vals)]
        if vals.size > 0:
            all_values_list.append(vals)

if len(all_values_list) == 0:
    raise RuntimeError("No valid numerical values ​​are available for plotting.")

all_values = np.concatenate(all_values_list)
y_min = 0.05
y_max_data = float(np.nanmax(all_values)) + 0.01
ytick_interval = 0.05 if (y_max_data - y_min) < 0.5 else 0.10
y_top_data = np.ceil(y_max_data / ytick_interval) * ytick_interval
y_range = y_top_data - y_min

axis_ylim_top = y_top_data + 0.08 * y_range

# figure
group_count   = len(dfs)
width_per_grp = 4.0
fig_width     = group_count * width_per_grp
fig_height    = fig_width * 0.40
plt.figure(figsize=(fig_width, fig_height))
ax = plt.gca()

# color
colors = [
    (21/255, 75/255, 168/255),    # G-G
    (92/255, 140/255, 199/255),   # G-W
    (1, 1, 1),                    # G-W/G-G
    (186/255, 205/255, 232/255),  # W-W
    (1, 1, 1)                     # W-W/G-G
]
border_colors = [
    (0, 0, 0),                    # G-G
    (21/255, 75/255, 168/255),    # G-W
    (21/255, 75/255, 168/255),    # G-W/G-G
    (92/255, 140/255, 199/255),   # W-W
    (92/255, 140/255, 199/255)    # W-W/G-G
]
correlation_labels = cond_order

within_group_spacing = 1.5   # Spacing between conditions within the same group
group_spacing        = 3.2   # Spacing between groups
pair_gap             = 0.20  
violin_width         = 0.60  
box_width            = 0.38  
cloud_alpha          = 0.8    

def draw_half_violin(ax, values, xcenter, face, edge, width=0.6, alpha=0.35):
    parts = ax.violinplot(values, positions=[xcenter], widths=width,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(face)
        pc.set_edgecolor(edge)
        pc.set_alpha(alpha)

        verts = pc.get_paths()[0].vertices
        verts[verts[:, 0] > xcenter, 0] = xcenter
        pc.set_zorder(1)
    return parts

center_positions = []
offset = 0.0
max_x = -np.inf

for i, df in enumerate(dfs):
    ds_name = dataset_order[i]
    group_centers = []
    for j, col in enumerate(df.columns):
        cond = col
        pos_center = offset + j * within_group_spacing
        group_centers.append(pos_center)

        vals = pd.to_numeric(df[col], errors='coerce').dropna().values
        if vals.size == 0:
            continue

        pos_cloud = pos_center - pair_gap
        pos_box   = pos_center + pair_gap

  
        draw_half_violin(ax, vals, pos_cloud,
                         face=colors[j], edge=border_colors[j],
                         width=violin_width, alpha=cloud_alpha)

  
        plt.boxplot(vals,
                    positions=[pos_box],
                    widths=box_width,
                    patch_artist=True,
                    boxprops=dict(facecolor=colors[j], edgecolor=border_colors[j], linewidth=1.6),
                    medianprops=dict(color=border_colors[j], linewidth=1.6),
                    whiskerprops=dict(color=border_colors[j], linewidth=1.4),
                    capprops=dict(color=border_colors[j], linewidth=1.4),
                    flierprops=dict(marker='o', markerfacecolor='none',
                                    markeredgecolor=border_colors[j], markersize=4, linestyle='none'))

       
        stars = stars_map.get(ds_name, {}).get(cond, "")
        if stars:
            local_max = float(np.nanmax(vals))
            y_star = local_max + 0.01 * y_range   
            ax.text(pos_center, y_star, stars,
                    ha='center', va='bottom',
                    fontsize=14, fontweight='bold')

        max_x = max(max_x, pos_box)

    center_positions.append(np.mean(group_centers))
    offset += within_group_spacing * len(df.columns) + group_spacing

# x-axis range
plt.xlim(-1, max_x + within_group_spacing)

# y-axis range & scale
plt.ylim(y_min, axis_ylim_top)
tick_start = 0.10
plt.yticks(np.arange(tick_start, y_top_data + 1e-9, ytick_interval), fontsize=12)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xticks(center_positions, name_list, ha='center', fontsize=18)
plt.ylabel("Prediction Accuracy", fontsize=18)

# figure legend
legend_patches = [plt.Rectangle((0, 0), 1, 1,
                                fc=colors[i], ec=border_colors[i],
                                lw=1.5, label=correlation_labels[i]) for i in range(5)]
plt.legend(handles=legend_patches,
           loc='upper left',
           bbox_to_anchor=(1.02, 1.0),
           borderaxespad=0,
           title="Legend",
           frameon=False)

plt.tight_layout()
plt.savefig(os.path.join(figure_folder, filename1), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(figure_folder, filename2), bbox_inches='tight')
plt.show()
