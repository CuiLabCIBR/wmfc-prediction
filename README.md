# WhiteMatter_Function_Prediction

Data and codes for our paper **"White-matter functional connectivity uniquely predicts brain maturation, cognition, and psychopathology beyond gray matter"**.

## Data Availability

The original data for these analyses are available via the following repositories:
* **ABCD-BIDS Community Collection (ABCC):** [NIMH Data Archive (NDA)](https://nda.nih.gov/edit_collection.html?id=3165)
* **Lifespan Human Connectome Project in Development (HCP-D)**
* **Reproducible Brain Charts (RBC):** [https://github.com/ReproBrainChart/PNC_BIDS](https://github.com/ReproBrainChart/PNC_BIDS)
* **Chinese Color Nest Project (devCCNP):** Science Data Bank (developmental component)
* **EFNY:** Data are not yet available, as collection is still ongoing.

## Software and System Requirements

The system requirements and installation guide for each software can be found on its respective website.

### Function & structural MRI processing
* **fMRIPrep v24.0.0** ([https://fmriprep.org/en/24.0.0/](https://fmriprep.org/en/24.0.0/))
* **XCP-D v0.7.1** ([https://xcp-d.readthedocs.io/en/0.7.1/](https://xcp-d.readthedocs.io/en/0.7.1/))
* **OS:** Linux

### Postprocessing
* **Connectome Workbench v1.5.0** ([https://www.humanconnectome.org/software/connectome-workbench](https://www.humanconnectome.org/software/connectome-workbench))
* **Python 3.9.18** ([https://www.python.org/](https://www.python.org/))
* **R v4.5.1** ([https://www.r-project.org](https://www.r-project.org))
* **MATLAB R2019a** ([https://www.mathworks.com/](https://www.mathworks.com/))
* **OS:** Linux

## Data

* **`demography`**: The demography folder contains the subject information used in this study.
* **`Prediction`**: The Prediction folder contains the correlation (*r*) and partial correlation (partial *r*) values for the age, cognition, and psychopathology prediction analyses.
* **`Atlas`**: The Atlas folder contains the parcellation files used in this study.

## Plot

Scripts for:
* Plotting group-level GM–GM, GM–WM, and WM–WM functional connectivity (FC) matrices.
* Visualizing prediction accuracies for age, cognition, and psychopathology across datasets and FC types.
* Visualizing network-level feature weights.

## step_01_data_processing

Scripts for data preprocessing and for screening participants based on head motion.

## step_02_calculate_FC

Scripts for computing GM–GM, GM–WM, and WM–WM FC and constructing the corresponding prediction features.

* `S01_segment_GM_WM.py`: Generates individual GM and WM tissue masks for each participant.
* `S02_resample_atlas.sh`: Resamples GM and WM atlases into each dataset's functional space.
* `S03_get_individualFCmatrix_*`: Extracts regional GM and WM time series and computes subject-level GM–GM, GM–WM, and WM–WM FC matrices.
* `S04_getZtrans_withTotalFC_*`: Applies Fisher z-transformation to FC values and concatenates them into participant-level prediction feature vectors.

## step_03_prediction

Scripts for performing partial least squares regression (PLS-R) and calculating prediction performance (Pearson’s *r*, MAE) and partial correlations (partial *r*).

* `S01_get_demo_*`: Prepares prediction targets and demographic/covariate information for each dataset.
* `PLSr1_CZ_Random_RegressCovariates_*`: Implements the PLS-R model with repeated cross-validation (CV) while regressing out covariates, and returns fold-wise predictions and model weights.
* `S02_Prediction_OverallPsyFactor_RandomCV_*`: Runs 5F-CV PLS-R prediction and permutation tests for age, cognition, and psychopathology using GM–GM, GM–WM, and WM–WM FC features, and saves prediction performance (Pearson’s *r*, MAE).
* `S03_get_totalPartialCorr_*`: Computes partial correlations between observed and predicted scores across FC types.

## step_04_get_significance

Scripts for evaluating the statistical significance of prediction performance and partial *r* using permutation tests.

* `S01_get_medianPartialCorr_perm_*`: Generates permutation-based null distributions of *r* and partial *r* across repetitions and FC types.
* `S02_get_medianPartialCorr_sig_*`: Computes *p*-values and significance for *r* and partial *r* and summarizes the results for visualization.

## step_05_feature_weights

Scripts for computing and visualizing Haufe-transformed feature weights, assessing their network-level significance, and generating Figures 2C, 3C, and 4B.

* `S01_getMedianWeight_*`: Computes median Haufe-transformed feature weights across CV repetitions for W–W, G–W, and G–G FC, separately for each prediction target and dataset, and generates the weight maps for Figures 3C and 4B.
* `S02_datasetsAveraged_*`: Averages median Haufe feature weights across datasets within each FC type to obtain cross-dataset weight patterns, and visualizes these patterns for Figure 2C.
* `S03_getSigNet_*`: Averages edge-wise Haufe weights into network-level scores and identifies significant network-level contributions.