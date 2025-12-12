%%
clear;
ProjectFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
targetStr_total = {'HCPD'; 'PNC'; 'EFNY'; 'CCNP'};

targetDir_total = {'/HCPD/code/4th_prediction/s02_PLSprediction/nosmooth/age_withHaufe_Schaefer100/HCPD_interview_age',...
    '/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100/PNC_age',...
    '/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth/age_withHaufe_Schaefer100_522/brainproject522_age',...
    '/CCNP/code/4th_prediction/nosmooth_n193_delete_under6/age_withHaufe_Schaefer100/CCNP_age'};

%% Initialize cell arrays to store results across different target variables
num_targets = length(targetStr_total);
permR_gg_totalStr = cell(num_targets, 3); %% targetString, corr, MAE
permR_gg_totalStr(:, 1) = targetStr_total;

permR_gw_totalStr = cell(num_targets, 3);
permR_gw_totalStr(:, 1) = targetStr_total;
permR_ww_totalStr = cell(num_targets, 3);
permR_ww_totalStr(:, 1) = targetStr_total;

perm_partialR_gw_totalStr = cell(num_targets, 2); %% targetString, corr
perm_partialR_gw_totalStr(:, 1) = targetStr_total;

perm_partialR_ww_totalStr = cell(num_targets, 2);
perm_partialR_ww_totalStr(:, 1) = targetStr_total;

num_cv_runs = 1000; % Number of permutation CV repetitions (Time_0 to Time_999)
num_folds = 5;     % Number of folds in the cross-validation

%% --- Main Loop for Each Target Variable ---
for i_str = 1:num_targets
    targetStr = targetStr_total{i_str};
    disp(['Processing target: ' targetStr]);
    targetDir = targetDir_total{i_str};
    BaseFolder = [ProjectFolder, '/', targetDir , '/RegressCovariates_RandomCV_Permutation'];

    % Initialize arrays for storing results across CV runs
    Corr_Overall_Perm_GG = nan(1, num_cv_runs);
    MAE_Overall_Perm_GG = nan(1, num_cv_runs);
    Corr_Overall_Perm_GW = nan(1, num_cv_runs);
    MAE_Overall_Perm_GW = nan(1, num_cv_runs);
    Corr_Overall_Perm_WW = nan(1, num_cv_runs);
    MAE_Overall_Perm_WW = nan(1, num_cv_runs);

   %% --- Load Overall Results (Corr, MAE) from Res_NFold.mat ---
   %% Gray matter-Gray matter (GGFC)
    fc_current = 'GGFC';
    disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
        if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Perm_GG(i) = tmp.Mean_Corr;
            MAE_Overall_Perm_GG(i) = tmp.Mean_MAE;
        else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
        end
    end
    permR_gg_totalStr{i_str, 2} = Corr_Overall_Perm_GG;
    permR_gg_totalStr{i_str, 3} = MAE_Overall_Perm_GG;

   %% Gray matter-White matter (GWFC)
    fc_current = 'GWFC';
    disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
         if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Perm_GW(i) = tmp.Mean_Corr;
            MAE_Overall_Perm_GW(i) = tmp.Mean_MAE;
         else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
         end
    end
    permR_gw_totalStr{i_str, 2} = Corr_Overall_Perm_GW;
    permR_gw_totalStr{i_str, 3} = MAE_Overall_Perm_GW;
  
   %% White matter-White matter (WWFC)
    fc_current = 'WWFC';
    disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
         if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Perm_WW(i) = tmp.Mean_Corr;
            MAE_Overall_Perm_WW(i) = tmp.Mean_MAE;
         else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
         end
    end
    permR_ww_totalStr{i_str, 2} = Corr_Overall_Perm_WW;
    permR_ww_totalStr{i_str, 3} = MAE_Overall_Perm_WW;
    
   %% --- Calculate Partial Correlations ---
    disp('  Calculating partial correlations...');
    perm_partialR_gw_total = nan(1, num_cv_runs);
    perm_partialR_ww_total = nan(1, num_cv_runs);

    for i_cv = 1:num_cv_runs
        i_cv
        % Get the *original* run index (0 to 100) for the i_cv-th best run
        % Handle potential NaN values in correlations if files were missing
        %%% sorted by CV times
        id_gg_current = i_cv; % This gives the index (1 to 101) in the original array
        id_gw_current = i_cv;
        id_ww_current = i_cv;

        % --- Load and Combine Scores from all 5 Folds for the current CV run ---
       %% GGFC Scores
        GG_index_allfolds = [];
        GG_predictScore_allfolds = [];
        GG_testScore_allfolds = [];
        valid_gg_run = true;
        for k = 0:(num_folds - 1) % Loop Folds 0 to 4
            FoldFile = fullfile(BaseFolder, ['Time_' num2str(id_gg_current - 1)], 'GGFC', ['Fold_' num2str(k) '_Score.mat']);
            if isfile(FoldFile)
                tmpFold = load(FoldFile);
                GG_index_allfolds = [GG_index_allfolds, tmpFold.Index];
                GG_predictScore_allfolds = [GG_predictScore_allfolds, tmpFold.Predict_Score];
                GG_testScore_allfolds = [GG_testScore_allfolds, tmpFold.Test_Score];
            else
                warning('File not found: %s. Partial correlation for rank %d might be inaccurate.', FoldFile, i_cv);
                valid_gg_run = false; break; % Stop loading folds for this run/type
            end
        end
        if ~valid_gg_run || isempty(GG_index_allfolds)
             warning('Could not load all GGFC fold files for run index %d (rank %d). Skipping partial correlation.', id_gg_current - 1, i_cv);
             continue; % Skip to next i_cv
        end
        [~, sort_idx_gg] = sort(GG_index_allfolds);
        GG_predictScore_allfolds_sorted = GG_predictScore_allfolds(sort_idx_gg);
        GG_testScore_allfolds_sorted = GG_testScore_allfolds(sort_idx_gg); % Keep original test scores aligned

       %% GWFC Scores
        GW_index_allfolds = [];
        GW_predictScore_allfolds = [];
        GW_testScore_allfolds = [];
        valid_gw_run = true;
        for k = 0:(num_folds - 1) % Loop Folds 0 to 4
            FoldFile = fullfile(BaseFolder, ['Time_' num2str(id_gw_current - 1)], 'GWFC', ['Fold_' num2str(k) '_Score.mat']);
             if isfile(FoldFile)
                tmpFold = load(FoldFile);
                GW_index_allfolds = [GW_index_allfolds, tmpFold.Index];
                GW_predictScore_allfolds = [GW_predictScore_allfolds, tmpFold.Predict_Score];
                GW_testScore_allfolds = [GW_testScore_allfolds, tmpFold.Test_Score];
             else
                warning('File not found: %s. Partial correlation for rank %d might be inaccurate.', FoldFile, i_cv);
                valid_gw_run = false; break;
             end
        end
         if ~valid_gw_run || isempty(GW_index_allfolds)
             warning('Could not load all GWFC fold files for run index %d (rank %d). Skipping partial correlation.', id_gw_current - 1, i_cv);
             continue; % Skip to next i_cv
         end
        [~, sort_idx_gw] = sort(GW_index_allfolds);
        GW_predictScore_allfolds_sorted = GW_predictScore_allfolds(sort_idx_gw);
        GW_testScore_allfolds_sorted = double(GW_testScore_allfolds(sort_idx_gw)); % Ensure double for corr

       %% WWFC Scores
        WW_index_allfolds = [];
        WW_predictScore_allfolds = [];
        WW_testScore_allfolds = [];
        valid_ww_run = true;
        for k = 0:(num_folds - 1) % Loop Folds 0 to 4
            FoldFile = fullfile(BaseFolder, ['Time_' num2str(id_ww_current - 1)], 'WWFC', ['Fold_' num2str(k) '_Score.mat']);
             if isfile(FoldFile)
                tmpFold = load(FoldFile);
                WW_index_allfolds = [WW_index_allfolds, tmpFold.Index];
                WW_predictScore_allfolds = [WW_predictScore_allfolds, tmpFold.Predict_Score];
                WW_testScore_allfolds = [WW_testScore_allfolds, tmpFold.Test_Score];
             else
                warning('File not found: %s. Partial correlation for rank %d might be inaccurate.', FoldFile, i_cv);
                valid_ww_run = false; break;
             end
        end
         if ~valid_ww_run || isempty(WW_index_allfolds)
             warning('Could not load all WWFC fold files for run index %d (rank %d). Skipping partial correlation.', id_ww_current - 1, i_cv);
             continue; % Skip to next i_cv
         end
        [~, sort_idx_ww] = sort(WW_index_allfolds);
        WW_predictScore_allfolds_sorted = WW_predictScore_allfolds(sort_idx_ww);
        WW_testScore_allfolds_sorted = double(WW_testScore_allfolds(sort_idx_ww)); % Ensure double for corr

        % --- Check if indices match across FC types (important for partial corr) ---
        if ~isequal(sort(GG_index_allfolds), sort(GW_index_allfolds)) || ~isequal(sort(GG_index_allfolds), sort(WW_index_allfolds))
            warning('Participant indices do not match across GG, GW, WW for rank %d run (IDs: %d, %d, %d). Skipping partial correlation.', i_cv, id_gg_current-1, id_gw_current-1, id_ww_current-1);
            continue; % Skip if participants aren't the same
        end

        % --- Perform Partial Correlation ---
        % Ensure scores are column vectors for partialcorr
        try
            [r_gw, ~] = partialcorr(GW_predictScore_allfolds_sorted(:), GW_testScore_allfolds_sorted(:), GG_predictScore_allfolds_sorted(:));
            [r_ww, ~] = partialcorr(WW_predictScore_allfolds_sorted(:), WW_testScore_allfolds_sorted(:), GG_predictScore_allfolds_sorted(:));
            perm_partialR_gw_total(1,i_cv) = r_gw;
            perm_partialR_ww_total(1,i_cv) = r_ww;
        catch ME
             warning('Error calculating partial correlation for rank %d (IDs: %d, %d, %d): %s. Skipping.', i_cv, id_gg_current-1, id_gw_current-1, id_ww_current-1, ME.message);
        end
    end
    % End loop over i_cv ranks
    perm_partialR_gw_totalStr{i_str, 2} = perm_partialR_gw_total;
    perm_partialR_ww_totalStr{i_str, 2} = perm_partialR_ww_total;

end
% End loop over target strings (i_str)

%% --- Save Results ---
disp('Saving results...');

% Prepare data cell for boxplot/further analysis (contains full distributions)
dataCell = cell(num_targets + 1, 6);
dataCell(1,:) = {'Target', 'PermCorr_GG', 'PermCorr_GW', 'PermPartialCorr_GW', 'PermCorr_WW', 'PermPartialCorr_WW'};
dataCell(2:end, 1) = targetStr_total;
dataCell(2:end, 2) = permR_gg_totalStr(:, 2);   % GG Correlations
dataCell(2:end, 3) = permR_gw_totalStr(:, 2);   % GW Correlations
dataCell(2:end, 4) = perm_partialR_gw_totalStr(:, 2); % GW Partial Correlations
dataCell(2:end, 5) = permR_ww_totalStr(:, 2);   % WW Correlations
dataCell(2:end, 6) = perm_partialR_ww_totalStr(:, 2); % WW Partial Correlations

% Define output file names (saving in the base project folder)
saveFolder = [ProjectFolder '/post_prediction/get_significance/Age'];
outputFileTotal = fullfile(saveFolder, 'perm_partial_results_total_5fold_sortedByCVtimes.mat');
outputFileBoxplot = fullfile(saveFolder, 'perm_partial_results_forBoxplot_5fold_sortedByCVtimes.mat');

% Save the comprehensive results
save(outputFileTotal, 'permR_gg_totalStr', 'permR_gw_totalStr','perm_partialR_gw_totalStr', ...
     'permR_ww_totalStr', 'perm_partialR_ww_totalStr', 'targetStr_total', '-v7.3'); % Use -v7.3 if variables might be large

% Save the data formatted for potential boxplots
save(outputFileBoxplot, 'dataCell', 'targetStr_total', '-v7.3');

disp(['Results saved to:']);
disp(['  ' outputFileTotal]);
disp(['  ' outputFileBoxplot]);
disp('Done.');
