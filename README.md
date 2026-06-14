# WhiteMatter_Function_Prediction

Data and codes for our paper **"White-matter functional connectivity uniquely predicts brain age, cognition, and psychopathology beyond gray matter"**.

## Data Availability

The original data for these analyses are available via the following repositories:
* **ABCD-BIDS Community Collection (ABCC):** [NIMH Data Archive (NDA)](https://nda.nih.gov/edit_collection.html?id=3165)
* **Lifespan Human Connectome Project in Development (HCP-D)**
* **Reproducible Brain Charts (RBC):** [https://github.com/ReproBrainChart/PNC_BIDS](https://github.com/ReproBrainChart/PNC_BIDS)
* **Chinese Color Nest Project (devCCNP):** Science Data Bank (developmental component)
* **EFNY:** Data are not yet available, as collection is still ongoing.

## Software and System Requirements

The system requirements and installation guide for each software can be found on its respective website.

### Functional and structural MRI processing
* **fMRIPrep v24.0.0** ([https://fmriprep.org/en/24.0.0/](https://fmriprep.org/en/24.0.0/))
* **XCP-D v0.7.1** ([https://xcp-d.readthedocs.io/en/0.7.1/](https://xcp-d.readthedocs.io/en/0.7.1/))
* **HCP Pipelines** ([https://github.com/Washington-University/HCPpipelines](https://github.com/Washington-University/HCPpipelines)) for the alternative EFNY workflow
* **Singularity/Apptainer** for running the fMRIPrep and XCP-D containers
* **Slurm** for the provided batch-submission scripts
* **OS:** Linux

### Postprocessing
* **Connectome Workbench v1.5.0** ([https://www.humanconnectome.org/software/connectome-workbench](https://www.humanconnectome.org/software/connectome-workbench))
* **Python 3.9.18** ([https://www.python.org/](https://www.python.org/))
* **R v4.5.1** ([https://www.r-project.org](https://www.r-project.org))
* **MATLAB R2019a** ([https://www.mathworks.com/](https://www.mathworks.com/))
* **OS:** Linux

## Data

* **`Data/demography/`**: Subject-level demographic, behavioral, and covariate tables used in the prediction analyses.
* **`Data/Prediction/`**: Summary prediction results used for reporting and visualization.
* **`Data/Atlas/`**: Parcellation images used to define the GM and WM regions.

## `step_01_data_processing/`

Scripts for MRI preprocessing, nuisance-regressor preparation, and participant/run screening:

* `S01_fmriprep_*.sh`: Dataset-specific Slurm jobs for running fMRIPrep.
* `S02_extract_confounds.py`: Extracts CSF and global-signal nuisance regressors from fMRIPrep confound tables.
* `S02_extract_confounds_by_title_ABCD.m`: Constructs and extracts ABCD nuisance regressors, including motion, tissue, and global signals.
* `S03_extract_columns_HCPD.m`: Standardizes and extracts the HCP-D nuisance-regressor columns required by XCP-D.
* `S03_xcpd_*.sh`: Dataset-specific Slurm jobs for running XCP-D with nuisance regression, temporal filtering, despiking, and no spatial smoothing.
* `S04_select_*`: Summarize framewise displacement and retain participants/runs satisfying the study's motion and time-point criteria.

### `step_01_data_processing/hcp_pipeline/`

Alternative EFNY preprocessing workflow based on the HCP Pipelines:

* `prepare_hcp_studyfolder_efny.py`: Converts EFNY BIDS inputs into the HCP-style StudyFolder layout using symbolic links.
* `PreFreeSurferPipelineBatch.sh`, `FreeSurferPipelineBatch.sh`, `PostFreeSurferPipelineBatch.sh`, `GenericfMRIVolumeProcessingPipelineBatch.sh`, and `GenericfMRISurfaceProcessingPipelineBatch.sh`: Wrappers for the structural and functional HCP Pipeline stages.
* `hcp_efny_env.sh` and `hcp_efny_batch_common.sh`: Shared environment configuration and stage helper functions.
* `submit_hcp_efny_chain.sh` and `submit_hcp_efny_stage.slurm.sh`: Submit ordered HCP Pipeline stages as dependent Slurm jobs.

## `step_02_calculate_FC/`

Scripts for computing GM–GM, GM–WM, and WM–WM FC and constructing the corresponding prediction features:

* `S01_segment_GM_WM.py`: Generates individual GM and WM tissue masks for each participant.
* `S02_resample_atlas.sh`: Resamples GM and WM atlases into each dataset's functional space.
* `S03_get_individualFCmatrix_*`: Extracts regional GM and WM time series and computes subject-level GM–GM, GM–WM, and WM–WM FC matrices.
* `S03_sbatch_individual_FC_*.sh`: Dataset-specific Slurm submission scripts for subject-level FC calculation.
* `S04_getZtrans_withTotalFC_*`: Applies Fisher z-transformation to FC values and concatenates them into participant-level prediction feature vectors.

## `step_03_prediction/`

Scripts for performing partial least squares regression (PLS-R) and calculating prediction performance (Pearson’s *r*, MAE) and partial correlations (partial *r*).

* `S01_get_demo_*`: Prepares prediction targets and demographic/covariate information for each dataset.
* `PLSr1_CZ_Random_RegressCovariates_*`: Implements the PLS-R model with repeated cross-validation (CV) while regressing out covariates, and returns fold-wise predictions and model weights.
* `S02_Prediction_OverallPsyFactor_RandomCV_*`: Runs 5F-CV PLS-R prediction and permutation tests for age, cognition, and psychopathology using GM–GM, GM–WM, and WM–WM FC features, and saves prediction performance (Pearson’s *r*, MAE).
* `S03_get_totalPartialCorr_*`: Computes partial correlations between observed and predicted scores across FC types.

### Validation and sensitivity analyses

* `V_feature_merge/`: Tests whether concatenating two or three FC feature classes improves prediction.
* `V_holdout/`: Runs half-sample holdout validation for age, cognition, and psychopathology. The ABCD analyses can use family-aware split definitions to prevent relatives from being divided across training and test samples.
* `V_siblings/`: Repeats ABCD cognition and psychopathology prediction after selecting one participant per family.
* `V_hcppipeline/`: Repeats EFNY age prediction using FC features generated from the alternative HCP Pipeline preprocessing workflow.

Each validation directory contains its own prediction entry scripts and, where needed, a modified PLS-R implementation.

## `step_04_get_significance/`

Scripts for evaluating the statistical significance of prediction performance and partial *r* using permutation tests.

* `S01_get_medianPartialCorr_perm_*`: Generates permutation-based null distributions of *r* and partial *r* across repetitions and FC types.
* `S02_get_medianPartialCorr_sig_*`: Computes *p*-values and significance for *r* and partial *r* and summarizes the results for visualization.

The `age`, `cognition`, and `pfactor` filename suffixes identify the prediction target analyzed by each script.

## `step_05_feature_weights/`

Scripts for computing and visualizing Haufe-transformed feature weights, assessing their network-level significance, and generating Figures 2C, 3C, and 4B.

* `S01a_getMedianWeight_ww_*`: Computes median Haufe-transformed WM–WM feature weights across CV repetitions.
* `S01b_getMedianWeight_gw_*`: Computes median Haufe-transformed GM–WM feature weights across CV repetitions.
* `S01c_getMedianWeight_gg_*`: Computes median Haufe-transformed GM–GM feature weights across CV repetitions.
* `S02a_ww_datasetsAveraged_age.m`, `S02b_gw_datasetsAveraged_age.m`, and `S02c_gg_datasetsAveraged_age.m`: Average age-prediction Haufe weights across datasets for each FC class.
* `S03_getSigNet_*`: Averages edge-wise Haufe weights into network-level scores and identifies significant network-level contributions.
* `JHU68_info.mat`, `Schaefer100_info.mat`, `gm_labels.mat`, and `wm_labels.mat`: Atlas labels, network assignments, and ordering information used by the feature-weight scripts.
* `cmap0.mat`: Diverging color map used for FC and feature-weight visualizations.

## `plot/`

Scripts for generating the main FC, prediction-performance, and feature-weight figures.
