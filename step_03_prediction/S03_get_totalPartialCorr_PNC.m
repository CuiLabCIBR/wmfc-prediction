%%
clear;
% targetStr_total = {'nihtbx_picvocab','nihtbx_picture_uncorrected','nihtbx_flanker','nihtbx_list', 'nihtbx_cardsort','nihtbx_pattern','nihtbx_reading'};

targetStr_total = {'PNC_age'};
targetStr_total = targetStr_total';

% Initialize cell arrays to store results across different target variables
num_targets = length(targetStr_total);
R_gg_totalStr = cell(num_targets, 3);
R_gg_totalStr(:, 1) = targetStr_total;
R_gw_totalStr = cell(num_targets, 3);
R_gw_totalStr(:, 1) = targetStr_total;
R_ww_totalStr = cell(num_targets, 3);
R_ww_totalStr(:, 1) = targetStr_total;

partialR_gw_totalStr = cell(num_targets, 2);
partialR_gw_totalStr(:, 1) = targetStr_total;
partialR_ww_totalStr = cell(num_targets, 2);
partialR_ww_totalStr(:, 1) = targetStr_total;

medianResults_totalStr = cell(num_targets, 6); % string,gg,gw,ww,gw_partial,ww_partial
medianResults_totalStr(:, 1) = targetStr_total;

% --- Define Project Folder ---
% !! Adjust this path if your base 'prediction_new/brain_project' is elsewhere !!
ProjectFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100';

num_cv_runs = 101; % Number of random CV repetitions (Time_0 to Time_100)
num_folds = 5;     % Number of folds in the cross-validation

%% --- Main Loop for Each Target Variable ---
for i_str = 1:num_targets
    targetStr = targetStr_total{i_str};
    disp(['Processing target: ' targetStr]);

    % Define the base directory for this target based on the new structure
    BaseFolder = fullfile(ProjectFolder, targetStr, 'RegressCovariates_RandomCV');
    if ~isfolder(BaseFolder)
        error('Base folder not found: %s', BaseFolder);
    end

    % Initialize arrays for storing results across CV runs
    Corr_Overall_Actual_GG = nan(1, num_cv_runs);
    MAE_Overall_Actual_GG = nan(1, num_cv_runs);
    Corr_Overall_Actual_GW = nan(1, num_cv_runs);
    MAE_Overall_Actual_GW = nan(1, num_cv_runs);
    Corr_Overall_Actual_WW = nan(1, num_cv_runs);
    MAE_Overall_Actual_WW = nan(1, num_cv_runs);

    %% --- Load Overall Results (Corr, MAE) from Res_NFold.mat ---
    % Gray matter-Gray matter (GGFC)
    fc_current = 'GGFC';
    disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
        if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Actual_GG(i) = tmp.Mean_Corr;
            MAE_Overall_Actual_GG(i) = tmp.Mean_MAE;
        else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
        end
    end
    % Get median results for GGFC
    Corr_median_GG = median(Corr_Overall_Actual_GG, 'omitnan');
    MAE_median_GG = median(MAE_Overall_Actual_GG, 'omitnan');
    disp(['    ' fc_current ' median corr: ' num2str(Corr_median_GG)]);
    medianResults_totalStr{i_str, 2} = Corr_median_GG;
    R_gg_totalStr{i_str, 2} = Corr_Overall_Actual_GG;
    R_gg_totalStr{i_str, 3} = MAE_Overall_Actual_GG;
    [~, ind_gg] = sort(Corr_Overall_Actual_GG, 'descend', 'MissingPlacement', 'last'); % Sort IDs by correlation
    median_id_gg = ind_gg(51)-1 % python starts with 0!!!

    % Gray matter-White matter (GWFC)
    fc_current = 'GWFC';
     disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
         if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Actual_GW(i) = tmp.Mean_Corr;
            MAE_Overall_Actual_GW(i) = tmp.Mean_MAE;
         else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
         end
    end
    % Get median results for GWFC
    Corr_median_GW = median(Corr_Overall_Actual_GW, 'omitnan');
    MAE_median_GW = median(MAE_Overall_Actual_GW, 'omitnan');
    disp(['    ' fc_current ' median corr: ' num2str(Corr_median_GW)]);
    medianResults_totalStr{i_str, 3} = Corr_median_GW;
    R_gw_totalStr{i_str, 2} = Corr_Overall_Actual_GW;
    R_gw_totalStr{i_str, 3} = MAE_Overall_Actual_GW;
    [~, ind_gw] = sort(Corr_Overall_Actual_GW, 'descend', 'MissingPlacement', 'last'); % Sort IDs by correlation
    median_id_gw = ind_gw(51)-1 % python starts with 0!!!
    
    % White matter-White matter (WWFC)
    fc_current = 'WWFC';
     disp(['  Loading overall results for ' fc_current '...']);
    for i = 1:num_cv_runs
        ResFile = fullfile(BaseFolder, ['Time_' num2str(i - 1)], fc_current, 'Res_NFold.mat');
         if isfile(ResFile)
            tmp = load(ResFile);
            Corr_Overall_Actual_WW(i) = tmp.Mean_Corr;
            MAE_Overall_Actual_WW(i) = tmp.Mean_MAE;
         else
            warning('File not found: %s. Skipping run %d for %s.', ResFile, i-1, fc_current);
         end
    end
    % Get median results for WWFC
    Corr_median_WW = median(Corr_Overall_Actual_WW, 'omitnan');
    MAE_median_WW = median(MAE_Overall_Actual_WW, 'omitnan');
    disp(['    ' fc_current ' median corr: ' num2str(Corr_median_WW)]);
    medianResults_totalStr{i_str, 4} = Corr_median_WW; % Note: Original code put this in col 4
    R_ww_totalStr{i_str, 2} = Corr_Overall_Actual_WW;
    R_ww_totalStr{i_str, 3} = MAE_Overall_Actual_WW;
    [~, ind_ww] = sort(Corr_Overall_Actual_WW, 'descend', 'MissingPlacement', 'last'); % Sort IDs by correlation
    median_id_ww = ind_ww(51)-1 % python starts with 0!!!
    
    %% --- Calculate Partial Correlations ---
    disp('  Calculating partial correlations...');
    partialR_gw_total = nan(1, num_cv_runs);
    partialR_ww_total = nan(1, num_cv_runs);

    for i_cv = 1:num_cv_runs
        % Get the *original* run index (0 to 100) for the i_cv-th best run
        % Handle potential NaN values in correlations if files were missing
        if i_cv > sum(~isnan(Corr_Overall_Actual_GG)) || i_cv > sum(~isnan(Corr_Overall_Actual_GW)) || i_cv > sum(~isnan(Corr_Overall_Actual_WW))
           warning('Skipping partial correlation calculation for rank %d due to missing data in overall results.', i_cv);
           continue; % Skip if this rank corresponds to a NaN correlation
        end
        %%% sorted by CORR 
%         id_gg_current = ind_gg(i_cv); % This gives the index (1 to 101) in the original array
%         id_gw_current = ind_gw(i_cv);
%         id_ww_current = ind_ww(i_cv);
        %%% sorted by CV times
        id_gg_current = i_cv; % This gives the index (1 to 101) in the original array
        id_gw_current = i_cv;
        id_ww_current = i_cv;

        % --- Load and Combine Scores from all 5 Folds for the current CV run ---

        % GGFC Scores
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

        % GWFC Scores
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

        % WWFC Scores
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
            partialR_gw_total(1, i_cv) = r_gw;
            partialR_ww_total(1, i_cv) = r_ww;
        catch ME
             warning('Error calculating partial correlation for rank %d (IDs: %d, %d, %d): %s. Skipping.', i_cv, id_gg_current-1, id_gw_current-1, id_ww_current-1, ME.message);
        end
    end % End loop over i_cv ranks

    % Store median partial correlations and the full distribution
    medianResults_totalStr{i_str, 5} = median(partialR_gw_total, 'omitnan');
    medianResults_totalStr{i_str, 6} = median(partialR_ww_total, 'omitnan');
    partialR_gw_totalStr{i_str, 2} = partialR_gw_total;
    partialR_ww_totalStr{i_str, 2} = partialR_ww_total;
    disp(['    Median partial corr GW (controlling for GG): ', num2str(medianResults_totalStr{i_str, 5})]);
    disp(['    Median partial corr WW (controlling for GG): ', num2str(medianResults_totalStr{i_str, 6})]);


end % End loop over target strings (i_str)

%% --- Save Results ---
disp('Saving results...');

% Create table for median results
medianResults_totalTable = cell2table(medianResults_totalStr);
colNames = {'Target', 'MedianCorr_GG', 'MedianCorr_GW', 'MedianCorr_WW', 'MedianPartialCorr_GW', 'MedianPartialCorr_WW'}; % Adjusted column names for clarity
medianResults_totalTable.Properties.VariableNames = colNames;
disp(medianResults_totalTable); % Display the median results table

% Prepare data cell for boxplot/further analysis (contains full distributions)
dataCell = cell(num_targets + 1, 6);
dataCell(1,:) = {'Target', 'Corr_GG_All', 'Corr_GW_All', 'PartialCorr_GW_All', 'Corr_WW_All', 'PartialCorr_WW_All'};
dataCell(2:end, 1) = targetStr_total;
dataCell(2:end, 2) = R_gg_totalStr(:, 2);   % GG Correlations
dataCell(2:end, 3) = R_gw_totalStr(:, 2);   % GW Correlations
dataCell(2:end, 4) = partialR_gw_totalStr(:, 2); % GW Partial Correlations
dataCell(2:end, 5) = R_ww_totalStr(:, 2);   % WW Correlations
dataCell(2:end, 6) = partialR_ww_totalStr(:, 2); % WW Partial Correlations

% Define output file names (saving in the base project folder)

saveFolder = ProjectFolder;
outputFileTotal = fullfile(saveFolder, 'partial_results_total_5fold_sortedByCVtimes.mat');
outputFileBoxplot = fullfile(saveFolder, 'partial_results_forBoxplot_5fold_sortedByCVtimes.mat');
% outputFileTotal = fullfile(saveFolder, 'partial_results_total_5fold_sortedByCORR.mat');
% outputFileBoxplot = fullfile(saveFolder, 'partial_results_forBoxplot_5fold_sortedByCORR.mat');

% Save the comprehensive results
save(outputFileTotal, 'partialR_gw_totalStr', 'partialR_ww_totalStr', ...
     'R_gg_totalStr', 'R_gw_totalStr', 'R_ww_totalStr', ...
     'medianResults_totalTable', 'targetStr_total', '-v7.3'); % Use -v7.3 if variables might be large

% Save the data formatted for potential boxplots
save(outputFileBoxplot, 'dataCell', 'targetStr_total', '-v7.3');

disp(['Results saved to:']);
disp(['  ' outputFileTotal]);
disp(['  ' outputFileBoxplot]);
disp('Done.');