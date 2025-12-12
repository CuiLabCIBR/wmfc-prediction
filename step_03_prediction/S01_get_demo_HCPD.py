#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple

TXT_SUBLIST = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/code/3rd_calFC/HCPD_subList_FDpassed_n531.txt')
AGE_SEX_CSV = Path('/ibmgpfs/cuizaixu_lab/zhaoshaoling/MSC_data/HCPD/code_WMpost/step04_HCPDprediction/step02_getAtlasFeature/networkFC/newConfounds/demoSelection/HCPD_sublist_age_sex_n652.csv')
SITE_CSV    = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/demography/HCPD_sublist_site.csv')
MOTION_CSV  = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/demography/HCPD_MotionInfo_n633_FDpassed_with4run_new.csv')

OUT_COVS    = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/code/4th_prediction/s01_preparedata/HCPD_covariates_531.csv')
OUT_COMB    = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/code/4th_prediction/s01_preparedata/HCPD_demoCombined_531.csv')
OUT_SITEMAP = Path('/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/code/4th_prediction/s01_preparedata/HCPD_site_code_map.csv')

OUT_TXT     = TXT_SUBLIST.with_name(TXT_SUBLIST.stem + '_siteNotNull.txt')


def read_subjects_from_txt(txt_path: Path) -> pd.DataFrame:
    lines = [ln.strip() for ln in txt_path.read_text(encoding='utf-8').splitlines() if ln.strip()]
    def norm(s: str) -> str:
        s = s.strip()
        if s.startswith('sub-'):
            s = s[4:]
        return s
    subs = [norm(x) for x in lines]
    seen = set(); ordered = []
    for s in subs:
        if s not in seen:
            ordered.append(s); seen.add(s)
    return pd.DataFrame({'subID': ordered})

def _auto_read_csv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, engine='python', sep=None)

def load_age_sex(path: Path) -> pd.DataFrame:
    df = _auto_read_csv(path)
    low = {c.lower(): c for c in df.columns}
    sid_col = low.get('src_subject_id') or low.get('subid') or low.get('subject') or low.get('sublist') or list(df.columns)[0]
    age_col = low.get('interview_age') or 'interview_age'
    sex_col = low.get('sex') or 'sex'
    df = df.rename(columns={sid_col: 'subID', age_col: 'interview_age', sex_col: 'sex'})
    df['subID'] = df['subID'].astype(str).str.replace(r'^sub-', '', regex=True)
    df['sex'] = df['sex'].astype(str).str.strip().str.upper().map({'M': 1, 'F': 0})
    return df[['subID', 'interview_age', 'sex']].copy()

def load_site(path: Path) -> pd.DataFrame:
    df = _auto_read_csv(path)
    cols_lower = {c.lower(): c for c in df.columns}
    sub_col = None; site_col = None
    for key in ['sublist', 'src_subject_id', 'subid', 'subject']:
        if key in cols_lower:
            sub_col = cols_lower[key]; break
    for key in ['site', 'sites', 'scanner', 'location']:
        if key in cols_lower:
            site_col = cols_lower[key]; break
    if sub_col is None:
        for c in df.columns:
            if 'sub' in c.lower():
                sub_col = c; break
    if site_col is None:
        for c in df.columns:
            if 'site' in c.lower():
                site_col = c; break
    if (sub_col is None or site_col is None) and df.shape[1] == 2:
        sub_col = sub_col or df.columns[0]
        site_col = site_col or df.columns[1]
    if sub_col is None or site_col is None:
        raise ValueError(f"Could not infer sublist/site columns in {path}. Columns found: {list(df.columns)}")
    df = df.rename(columns={sub_col: 'subID', site_col: 'site'})
    df['subID'] = df['subID'].astype(str).str.replace(r'^sub-', '', regex=True)
    df['site'] = df['site'].astype('string').str.strip()
    return df[['subID', 'site']].copy()

def load_fd(path: Path) -> pd.DataFrame:
    df = _auto_read_csv(path)
    cols_lower = {c.lower(): c for c in df.columns}
    sid_col = cols_lower.get('sublist') or cols_lower.get('src_subject_id') or cols_lower.get('subid') or list(df.columns)[0]
    meanfd_col = None
    for key in ['meanfd', 'mean_fd']:
        if key in cols_lower:
            meanfd_col = cols_lower[key]; break
    if meanfd_col is None:
        for c in df.columns:
            if 'mean' in c.lower() and 'fd' in c.lower():
                meanfd_col = c; break
    if meanfd_col is None:
        raise ValueError(f"Could not find meanFD column in {path}. Columns found: {list(df.columns)}")
    df = df.rename(columns={sid_col: 'subID', meanfd_col: 'meanFD'})
    df['subID'] = df['subID'].astype(str).str.replace(r'^sub-', '', regex=True)
    return df[['subID', 'meanFD']].copy()

def make_site_codes(site_series: pd.Series, fixed_map: Dict[str, int] = None) -> Tuple[pd.Series, Dict[str, int]]:
    sites = site_series.dropna().astype(str).str.strip()
    if fixed_map:
        mapping = fixed_map.copy()
    else:
        unique_sites = sorted(sites.unique())
        mapping = {name: i+1 for i, name in enumerate(unique_sites)}
    coded = site_series.map(mapping)
    return coded, mapping

def main():
    for p in [TXT_SUBLIST, AGE_SEX_CSV, SITE_CSV, MOTION_CSV]:
        if not p.exists():
            sys.exit(f"[ERROR] Missing input file: {p}")

    df_ids = read_subjects_from_txt(TXT_SUBLIST)
    df_age = load_age_sex(AGE_SEX_CSV)
    df_site = load_site(SITE_CSV)
    df_fd = load_fd(MOTION_CSV)

    # Merge
    df = (df_ids.merge(df_age, on='subID', how='left')
                 .merge(df_site, on='subID', how='left')
                 .merge(df_fd, on='subID', how='left'))

    # === Filter subjects whose site is empty (NA/empty string/'nan') and export a new TXT list. ===
    site_str = df['site'].astype('string')
    mask_site_ok = site_str.notna() & (site_str.str.strip() != '') & (site_str.str.lower() != 'nan')
    n_removed = (~mask_site_ok).sum()
    df = df.loc[mask_site_ok].copy()

    # Save the filtered participant list in TXT format.
    OUT_TXT.write_text('\n'.join(df['subID'].tolist()) + '\n', encoding='utf-8')

    # Subsequent exports are all based on the filtered data.
    site_code_series, site_map = make_site_codes(df['site'], fixed_map=None)

    covs = df[['subID', 'sex', 'meanFD']].copy()
    covs['site'] = site_code_series.astype('Int64')
    covs.to_csv(OUT_COVS, index=False)

    demo = df[['subID', 'interview_age', 'sex', 'meanFD', 'site']].copy()
    demo.to_csv(OUT_COMB, index=False)

    pd.DataFrame({'site': list(site_map.keys()), 'site_code': list(site_map.values())}) \
        .sort_values('site_code') \
        .to_csv(OUT_SITEMAP, index=False)

    # reports
    n_total = len(df_ids)
    n_covs = len(covs)
    n_demo = len(demo)
    missing_age = demo['interview_age'].isna().sum()
    missing_sex = demo['sex'].isna().sum()
    missing_fd = demo['meanFD'].isna().sum()

    print(f"Subjects in TXT (input): {n_total}")
    print(f"Removed with empty site: {n_removed}")
    print(f"After filter -> rows: {len(df)} (covariates: {n_covs}; combined: {n_demo})")
    print(f"Missing counts (after filter) -> age: {missing_age}, sex: {missing_sex}, meanFD: {missing_fd}")
    print(f"Filtered subject list saved to: {OUT_TXT.resolve()}")
    print(f"Site code map saved to: {OUT_SITEMAP.resolve()}")
    print(f"Wrote:\n  - {OUT_COVS.resolve()}\n  - {OUT_COMB.resolve()}")

if __name__ == '__main__':
    main()