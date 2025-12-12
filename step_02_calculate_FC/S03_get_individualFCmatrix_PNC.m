function s02_get_individualFCmatrix(subj)
    rest_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results/smooth/save_nii_final';
    work_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/code/3rd_calFC';
    results_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results'
    GM_region_num = 100; 
    WM_region_num = 68;
    % WM_region_num = 50;

    %% make individual FC dir
    subj_path = [results_dir '/FCmatrix_individual_schaefer100_nosmooth/' subj];
    mkdir(subj_path);

    %% gray matter
    gm_bold_data = niftiread([rest_dir '/'  subj '/ses-PNC1/func/' subj '_ses-PNC1_task-rest_acq-singleband_space-MNI152NLin6Asym_res-2_desc-denoised_bold_GM.nii.gz']);
    %%%1. get GM atlas
    gm_atlas_data = niftiread([work_dir '/atlas/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm_resliced.nii']);
    gm_atlas_matrix = gm_atlas_data;
    %%% 2. GG FC compute with atlas
    gm_matrix_totalNetwork_averaged = [];           
    for i_net =1:GM_region_num
        index = find(gm_atlas_matrix==i_net);
        gm_matrix_currentNetwork = [];
        for i_timepoint = 1:size(gm_bold_data, 4)
            current_timepoint = gm_bold_data(:, :, :, i_timepoint);
            gm_matrix_currentNetwork (:, i_timepoint) = current_timepoint(index);
        end
        gm_matrix_currentNetwork_averaged = mean(gm_matrix_currentNetwork);
        gm_matrix_totalNetwork_averaged(i_net,:) = gm_matrix_currentNetwork_averaged;
    end
     GG_FCmatrix = corr(gm_matrix_totalNetwork_averaged');
     save([subj_path '/' subj '_Schaefer100_tiemseries_gm.mat'], 'gm_matrix_totalNetwork_averaged');
     save([subj_path '/' subj '_Schaefer100_GG_FC.mat'], 'GG_FCmatrix');
     clear gm_bold_data;
    
    %% white matter
    wm_bold_data = niftiread([rest_dir '/'  subj '/ses-PNC1/func/' subj '_ses-PNC1_task-rest_acq-singleband_space-MNI152NLin6Asym_res-2_desc-denoised_bold_WM.nii.gz']);
    %%%1. get GM atlas
    wm_atlas_data = niftiread([work_dir '/atlas/rICBM_DTI_81_WMPM_60p_FMRIB58_resliced.nii']);
    % wm_atlas_data = niftiread([work_dir '/atlas/rICBM_DTI_81_WMPM_FMRIB58_resliced_withniftiwrite.nii']);
%     wm_atlas_matrix = wm_atlas_data.img;
    wm_atlas_matrix = wm_atlas_data;
    %%% 2. WW FC compute with atlas
    wm_matrix_totalNetwork_averaged = [];           
    for i_net =1:WM_region_num
        index = find(wm_atlas_matrix==i_net);
        wm_matrix_currentNetwork = [];
        for i_timepoint = 1:size(wm_bold_data, 4)
            current_timepoint = wm_bold_data(:, :, :, i_timepoint);
            wm_matrix_currentNetwork (:, i_timepoint) = current_timepoint(index);
        end
        wm_matrix_currentNetwork_averaged = mean(wm_matrix_currentNetwork);
        wm_matrix_totalNetwork_averaged(i_net,:) = wm_matrix_currentNetwork_averaged;
    end
     WW_FCmatrix = corr(wm_matrix_totalNetwork_averaged');
     save([subj_path '/' subj '_rICBM_60p_tiemseries_wm.mat'], 'wm_matrix_totalNetwork_averaged');
     save([subj_path '/' subj '_rICBM_60p_WW_FC.mat'], 'WW_FCmatrix');
     % save([subj_path '/' subj '_rICBM_50_tiemseries_wm.mat'], 'wm_matrix_totalNetwork_averaged');
     % save([subj_path '/' subj '_rICBM_50_WW_FC.mat'], 'WW_FCmatrix');
     
    %% compute GW FC
    for i = 1:size(gm_matrix_totalNetwork_averaged, 1)
        for j = 1:size(wm_matrix_totalNetwork_averaged, 1)          
            GW_FCmatrix(i, j) = corr(gm_matrix_totalNetwork_averaged(i, :)', wm_matrix_totalNetwork_averaged(j, :)');
        end
    end
    save([subj_path '/' subj '_Schaefer100_rICBM_60p_GW_FC.mat'], 'GW_FCmatrix');

end
