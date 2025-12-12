function s02_get_HCPD_indiviFC_4run(subj)
    % clear;
    % subj = 'sub-HCD0001305'; 
        smooth_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/results/smooth';
        output_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/results';
        work_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/code/3rd_calFC';
        GM_region_num = 100; 
        WM_region_num = 68;
    
        %% make individual FC dir
        subj_path = [output_dir '/FCmatrix_individual_Schaefer100_nosmooth/' subj];
        mkdir(subj_path);
    
        %% gray matter
        gm_bold_data1 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST1_acq-AP_space-MNI152NLin6Asym_desc-denoised_bold_GM.nii.gz']);
        gm_bold_data2 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST1_acq-PA_space-MNI152NLin6Asym_desc-denoised_bold_GM.nii.gz']);
        gm_bold_data3 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST2_acq-AP_space-MNI152NLin6Asym_desc-denoised_bold_GM.nii.gz']);
        gm_bold_data4 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST2_acq-PA_space-MNI152NLin6Asym_desc-denoised_bold_GM.nii.gz']);
        gm_bold_data = cat(4, gm_bold_data1, gm_bold_data2, gm_bold_data3, gm_bold_data4);
        clear gm_bold_data1 gm_bold_data2 gm_bold_data3 gm_bold_data4
        %%%1. get GM atlas
        gm_atlas_data = niftiread([work_dir '/atlas/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm_resliced_new.nii']);
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
    %      save([subj_path '/' subj '_Yeo7_region51_tiemseries_gm_4runCombined.mat'], 'gm_matrix_totalNetwork_averaged');
    %      save([subj_path '/' subj '_Yeo7_region51_GG_FC.mat'], 'GG_FCmatrix');
         save([subj_path '/' subj '_Schaefer100_tiemseries_gm_4runCombined.mat'], 'gm_matrix_totalNetwork_averaged');
         save([subj_path '/' subj '_Schaefer100_GG_FC.mat'], 'GG_FCmatrix');
         clear gm_bold_data;
        
        %% white matter
        wm_bold_data1 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST1_acq-AP_space-MNI152NLin6Asym_desc-denoised_bold_WM.nii.gz']);
        wm_bold_data2 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST1_acq-PA_space-MNI152NLin6Asym_desc-denoised_bold_WM.nii.gz']);
        wm_bold_data3 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST2_acq-AP_space-MNI152NLin6Asym_desc-denoised_bold_WM.nii.gz']);
        wm_bold_data4 = niftiread([smooth_dir '/'  subj '/func/' subj '_task-REST2_acq-PA_space-MNI152NLin6Asym_desc-denoised_bold_WM.nii.gz']);
        wm_bold_data = cat(4, wm_bold_data1, wm_bold_data2, wm_bold_data3, wm_bold_data4);
        clear wm_bold_data1 wm_bold_data2 wm_bold_data3 wm_bold_data4
        %%%1. get WM atlas
        wm_atlas_data = niftiread([work_dir '/atlas/rICBM_DTI_81_WMPM_60p_FMRIB58_resliced_new.nii']);
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
    %      save([subj_path '/' subj '_rICBM_region50_tiemseries_wm_4runCombine.mat'], 'wm_matrix_totalNetwork_averaged');
    %      save([subj_path '/' subj '_rICBM_region50_WW_FC.mat'], 'WW_FCmatrix');
         save([subj_path '/' subj '_rICBM_60p_tiemseries_wm_4runCombine.mat'], 'wm_matrix_totalNetwork_averaged');
         save([subj_path '/' subj '_rICBM_60p_WW_FC.mat'], 'WW_FCmatrix');
         clear wm_bold_data;
         
        %% compute GW FC
        for i = 1:size(gm_matrix_totalNetwork_averaged, 1)
            for j = 1:size(wm_matrix_totalNetwork_averaged, 1)          
                GW_FCmatrix(i, j) = corr(gm_matrix_totalNetwork_averaged(i, :)', wm_matrix_totalNetwork_averaged(j, :)');
            end
        end
    %     save([subj_path '/' subj '_Yeo7_rICBM_GW_FC.mat'], 'GW_FCmatrix');
        save([subj_path '/' subj '_Schaefer100_rICBM_60p_GW_FC.mat'], 'GW_FCmatrix');
    
    end
    