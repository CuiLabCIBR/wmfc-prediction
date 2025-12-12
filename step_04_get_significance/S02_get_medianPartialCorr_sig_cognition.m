%%
clear;
ProjectFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
targetStr_total = {'nihtbx_totalcomp'; 'nihtbx_fluidcomp'; 'nihtbx_cryst'};

resultsFolder= [ProjectFolder '/post_prediction/get_significance/Cognition']; 
load([resultsFolder '/perm_partial_results_total_5fold_sortedByCVtimes.mat']);
load([ProjectFolder '/ABCD/code/5th_prediction/PLS_prediction/networkFC/cognition_nosmooth_motion2runFD/partial_results_forBoxplot_5fold_sortedByCVtimes.mat']);

%% get perm results for GG,GW,WW 
for i_str = 1:length(targetStr_total)
    thisTarget = targetStr_total{i_str};
   %% GG
    MedianCorr_GG = median(dataCell{i_str+1, 2});       %obs:The median of Mean_Corr in 101 actual prediction
    PermCorr_GG = permR_gg_totalStr{i_str, 2}; %% Mean_Corr,1 is target string, 3 is MAE
    nan_idx = find(~isfinite(PermCorr_GG));
    
    if ~isempty(nan_idx)
        error('[NaN DETECTED] Target=%s | GG perm NaN count=%d | TimeIDs(0-based)=%s', ...
              thisTarget, numel(nan_idx), mat2str(nan_idx(:)'-1));
    end
    if numel(PermCorr_GG) ~= 1000
        error('[LEN MISMATCH] Target=%s | GG perm length=%d (expected 1000)', ...
              thisTarget, numel(PermCorr_GG));
    end

    if MedianCorr_GG>0
        perm_P_GG = length(find(PermCorr_GG >= MedianCorr_GG)) / 1000;
    elseif MedianCorr_GG<0
        perm_P_GG = length(find(PermCorr_GG <= MedianCorr_GG)) / 1000;
    else
        error('MedianCorr_GG == 0, Target=%s', thisTarget);
    end
    perm_P_GG_totalStr(i_str,1) = perm_P_GG;
     
   %% GW
    MedianCorr_GW = median(dataCell{i_str+1, 3});
    PermCorr_GW = permR_gw_totalStr{i_str, 2}; %% 1 is target string, 3 is MAE

    nan_idx = find(~isfinite(PermCorr_GW));
    if ~isempty(nan_idx)
        error('[NaN DETECTED] Target=%s | GW perm NaN count=%d | TimeIDs(0-based)=%s', ...
              thisTarget, numel(nan_idx), mat2str(nan_idx(:)'-1));
    end
    if numel(PermCorr_GW) ~= 1000
        error('[LEN MISMATCH] Target=%s | GW perm length=%d (expected 1000)', ...
              thisTarget, numel(PermCorr_GW));
    end

    if MedianCorr_GW>0
        perm_P_GW = length(find(PermCorr_GW >= MedianCorr_GW)) / 1000;
    elseif MedianCorr_GW<0
        perm_P_GW = length(find(PermCorr_GW <= MedianCorr_GW)) / 1000;
    else
        error('MedianCorr_GW == 0, Target=%s', thisTarget);
    end
    perm_P_GW_totalStr(i_str,1) = perm_P_GW;
    
   %% WW
    MedianCorr_WW = median(dataCell{i_str+1, 5}); 
    PermCorr_WW = permR_ww_totalStr{i_str, 2}; %% 1 is target string, 3 is MAE

    nan_idx = find(~isfinite(PermCorr_WW));
    if ~isempty(nan_idx)
        error('[NaN DETECTED] Target=%s | WW perm NaN count=%d | TimeIDs(0-based)=%s', ...
              thisTarget, numel(nan_idx), mat2str(nan_idx(:)'-1));
    end
    if numel(PermCorr_WW) ~= 1000
       error('[LEN MISMATCH] Target=%s | WW perm length=%d (expected 1000)', ...
              thisTarget, numel(PermCorr_WW));
    end

    if MedianCorr_WW>0
        perm_P_WW = length(find(PermCorr_WW >= MedianCorr_WW)) / 1000;
    elseif MedianCorr_WW<0
        perm_P_WW = length(find(PermCorr_WW <= MedianCorr_WW)) / 1000;
    else
        error('MedianCorr_WW == 0, Target=%s', thisTarget);       
    end
    perm_P_WW_totalStr(i_str,1) = perm_P_WW;
    
end

%% get perm results for partialGW and partialWW 
for i_str = 1:length(targetStr_total)
    thisTarget = targetStr_total{i_str};
    %% partial GW | GG
    MedianPartialCorr_GW = median(dataCell{i_str+1, 4}); 
    PermPartialCorr_GW = perm_partialR_gw_totalStr{i_str, 2}; 

    nan_idx = find(~isfinite(PermPartialCorr_GW));
    if ~isempty(nan_idx)
        error('[NaN DETECTED] Target=%s | partial GW perm NaN count=%d | TimeIDs(0-based)=%s', ...
              thisTarget, numel(nan_idx), mat2str(nan_idx(:)'-1));
    end
    if numel(PermPartialCorr_GW) ~= 1000
        error('[LEN MISMATCH] Target=%s | partial GW perm length=%d (expected 1000)', ...
              thisTarget, numel(PermPartialCorr_GW));
    end

    if MedianPartialCorr_GW>0
        perm_partialP_GW = length(find(PermPartialCorr_GW >= MedianPartialCorr_GW)) / 1000;
    elseif MedianPartialCorr_GW<0
        perm_partialP_GW = length(find(PermPartialCorr_GW <= MedianPartialCorr_GW)) / 1000;
    else
        error('MedianPartialCorr_GW == 0, Target=%s', thisTarget);
    end
    perm_partialP_GW_totalStr(i_str,1) = perm_partialP_GW;

   %% partial WW | GG
    MedianPartialCorr_WW = median(dataCell{i_str+1, 6});
    PermPartialCorr_WW = perm_partialR_ww_totalStr{i_str, 2}; %% 1 is target string, 3 is MAE
    
    nan_idx = find(~isfinite(PermPartialCorr_WW));
    if ~isempty(nan_idx)
        error('[NaN DETECTED] Target=%s | partial WW perm NaN count=%d | TimeIDs(0-based)=%s', ...
              thisTarget, numel(nan_idx), mat2str(nan_idx(:)'-1));
    end
    if numel(PermPartialCorr_WW) ~= 1000
        error('[LEN MISMATCH] Target=%s | partial WW perm length=%d (expected 1000)', ...
              thisTarget, numel(PermPartialCorr_WW));
    end

    if MedianPartialCorr_WW>0
        perm_partialP_WW = length(find(PermPartialCorr_WW >= MedianPartialCorr_WW)) / 1000;
    elseif MedianPartialCorr_WW<0
        perm_partialP_WW = length(find(PermPartialCorr_WW <= MedianPartialCorr_WW)) / 1000;
    else
        error('MedianPartialCorr_WW == 0, Target=%s（不应出现）', thisTarget);
    end
    perm_partialP_WW_totalStr(i_str,1) = perm_partialP_WW;
    
end

%% merge p value
P_total = cell(length(targetStr_total)+1,6);
P_total(1,:) = {'Target', 'P_GG', 'P_GW', 'P_WW', 'partialP_GW', 'partialP_WW'};
for i_str = 1:length(targetStr_total)
    currentStr = targetStr_total{i_str};
    P_total{i_str+1,1} = currentStr;
    P_current = [perm_P_GG_totalStr(i_str), perm_P_GW_totalStr(i_str), perm_P_WW_totalStr(i_str)];
    partialP_current = [perm_partialP_GW_totalStr(i_str), perm_partialP_WW_totalStr(i_str)];
    P_total(i_str+1,2:6) = num2cell([P_current, partialP_current]);
end
   

%% FDR（Benjamini–Hochberg） across all datasets × targets（4×5=20）
p_mat = nan(length(targetStr_total), 5);
for i_str = 1:length(targetStr_total)
    p_mat(i_str, :) = [ ...
        perm_P_GG_totalStr(i_str), ...
        perm_P_GW_totalStr(i_str), ...
        perm_P_WW_totalStr(i_str), ...
        perm_partialP_GW_totalStr(i_str), ...
        perm_partialP_WW_totalStr(i_str) ];
end

% Flatten, filter out non-finite values and sort in ascending order.
p_vec  = p_mat(:);
valid  = isfinite(p_vec);
ps     = p_vec(valid);
[ps_sorted, order] = sort(ps);
m = numel(ps_sorted);

% BH step-up
q_sorted    = nan(size(ps_sorted));
q_sorted(m) = ps_sorted(m);
for k = m-1:-1:1
    q_sorted(k) = min(ps_sorted(k) * m / k, q_sorted(k+1));
end
q_sorted = min(q_sorted, 1);  

% according to the original position
q_unsorted        = nan(size(ps));
q_unsorted(order) = q_sorted;

% Re-fill the vector and restore the matrix shape
q_vec        = nan(size(p_vec));
q_vec(valid) = q_unsorted;
q_mat        = reshape(q_vec, size(p_mat));

%  Q_total
Q_total = cell(length(targetStr_total)+1, 6);
Q_total(1,:) = {'Target', 'q_GG', 'q_GW', 'q_WW', 'q_partial_GW', 'q_partial_WW'};
for i_str = 1:length(targetStr_total)
    Q_total{i_str+1,1}   = targetStr_total{i_str};
    Q_total(i_str+1,2:6) = num2cell(q_mat(i_str,:));
end

% save
save([resultsFolder '/Pvalue_byPermutation.mat'], 'P_total', 'Q_total');
save([resultsFolder '/Pvalue_all.mat'], ...
     'P_total', 'Q_total', ...
     'perm_P_GG_totalStr', 'perm_P_GW_totalStr', 'perm_P_WW_totalStr', ...
     'perm_partialP_GW_totalStr', 'perm_partialP_WW_totalStr');

disp('Saved:');
disp([resultsFolder '/Pvalue_byPermutation.mat']);
disp([resultsFolder '/Pvalue_all.mat']);


