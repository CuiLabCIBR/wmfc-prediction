%% Permutation test for block-averaged weights: W-W, G-W, G-G
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');

datasets_str = {'HCPD', 'PNC', 'EFNY', 'CCNP'};
datasets_dir = { ...
    [proj_dir '/HCPD/code/4th_prediction/s02_PLSprediction/nosmooth/age_withHaufe_Schaefer100/HCPD_interview_age'], ... % HCPD
    [proj_dir '/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100/PNC_age'], ... % PNC
    [proj_dir '/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth/age_withHaufe_Schaefer100_522/brainproject522_age'], ... % EFNY
    [proj_dir '/CCNP/code/4th_prediction/nosmooth_n193_delete_under6/age_withHaufe_Schaefer100/CCNP_age'] ... % CCNP
};

nSet  = length(datasets_str);
out_dir = [proj_dir '/Figure/Age_perdiction_feature_weight'];

% find permutation number
perm_dirs_all = cell(nSet, 1);
for i_set = 1:nSet
    perm_dirs_all{i_set} = g_ls([datasets_dir{i_set} '/RegressCovariates_RandomCV_Permutation/Time_*']);
end
nPerm = length(perm_dirs_all{1});
fprintf('Detected %d permutations.\n', nPerm);

for i_set = 2:nSet
    if length(perm_dirs_all{i_set}) ~= nPerm
        error('Permutation number mismatch between datasets (%s vs %s).', ...
            datasets_str{1}, datasets_str{i_set});
    end
end

%% ---------- mean real data（WW / GW / GG） ----------
load wm_labels.mat;   % 68x1, 1-5
load gm_labels.mat;   % 100x1, 1-7

% ---- W-W: 68x68 ----
realW_mat = load([out_dir '/datasetsAveraged_Haufe_FCmatrix_ww.mat'], 'FCmatrix_full');
Mat_WW = realW_mat.FCmatrix_full;

Net_real_WW = zeros(5, 5);
for w_row = 1:5
    for w_col = 1:5
        idx_r = (wm_labels == w_row);
        idx_c = (wm_labels == w_col);
        subM  = Mat_WW(idx_r, idx_c);
        Net_real_WW(w_row, w_col) = mean(subM(:));
    end
end

% ---- G-W: 100x68 ----
realGW_mat = load([out_dir '/datasetsAveraged_Haufe_FCmatrix_gw.mat'], 'FCmatrix_full');
Mat_GW = realGW_mat.FCmatrix_full;     % 100 x 68

Net_real_GW = zeros(7, 5);
for g_row = 1:7
    for w_col = 1:5
        idx_g = (gm_labels == g_row);  % 100
        idx_w = (wm_labels == w_col);  % 68
        subM  = Mat_GW(idx_g, idx_w);
        Net_real_GW(g_row, w_col) = mean(subM(:));
    end
end

% ---- G-G: 100x100 ----
realGG_mat = load([out_dir '/datasetsAveraged_Haufe_FCmatrix_gg.mat'], 'FCmatrix_full');
Mat_GG = realGG_mat.FCmatrix_full;

Net_real_GG = zeros(7, 7);
for g_row = 1:7
    for g_col = 1:7
        idx_r = (gm_labels == g_row);
        idx_c = (gm_labels == g_col);
        subM  = Mat_GG(idx_r, idx_c);
        Net_real_GG(g_row, g_col) = mean(subM(:));
    end
end

%% ---------- Part 1: W-W permutation（68×68 → 5×5） ----------
fprintf('Permutation for W-W blocks...\n');

Net_perm_WW = nan(5, 5, nPerm);  % 5x5xNperm

for i_perm = 1:nPerm
    Net_sets = nan(5, 5, nSet); 
    
    for i_set = 1:nSet
        current_set = datasets_str{i_set};
        data_dir    = datasets_dir{i_set};
        this_perm_dir = perm_dirs_all{i_set}{i_perm};

        % 1) Aggregate the Haufe weights of all folds under this permutation.
        fold_files = g_ls([this_perm_dir '/WWFC/Fold_*_Score.mat']);
        clear W_cv
        for i_cv = 1:length(fold_files)
            tmp = load(fold_files{i_cv});
            W_cv(i_cv, :) = tmp.w_Brain_Haufe(:)';  
        end
        w_perm = mean(W_cv, 1);   % 1×nEdges

        % 2) vector -> 68x68 
        lower_tri_vec = w_perm;
        m = 68; n = 68;
        if numel(lower_tri_vec) ~= m*(m-1)/2
            error('WW vector length mismatch in %s perm %d', current_set, i_perm);
        end

        FC = zeros(m, n);
        if strcmpi(current_set, 'EFNY')
            % Python 
            [ri, ci] = find(tril(true(m, n), -1));
            [~, ord] = sortrows([ri, ci], [1 2]);
            lin_idx   = sub2ind([m n], ri(ord), ci(ord));
            FC(lin_idx) = lower_tri_vec(:);
        else
            % MATLAB 
            tril_idx = find(tril(true(m, n), -1));
            FC(tril_idx) = lower_tri_vec(:);
        end
        FC_full = FC + FC';
        FC_full(logical(eye(m))) = 0;

        % 3) 68x68 → 5x5 net
        Net_tmp = zeros(5,5);
        for w_row = 1:5
            for w_col = 1:5
                idx_r = (wm_labels == w_row);
                idx_c = (wm_labels == w_col);
                subM  = FC_full(idx_r, idx_c);
                Net_tmp(w_row, w_col) = mean(subM(:));
            end
        end

        Net_sets(:,:,i_set) = Net_tmp;
    end

    % 4) average fou datasets
    Net_perm_WW(:,:,i_perm) = mean(Net_sets, 3);
end

%% cal W-W p value
p_WW     = nan(5,5);
pFDR_WW  = nan(5,5);
for i = 1:5
    for j = 1:5
        real_val  = Net_real_WW(i,j);
        null_vals = squeeze(Net_perm_WW(i,j,:));

        if real_val > 0
            p_WW(i,j) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_WW(i,j) = sum(null_vals <= real_val) / nPerm;
        else
            p_WW(i,j) = NaN;
        end
    end
end
valid_idx = ~isnan(p_WW(:));
p_tmp = nan(size(p_WW));
p_tmp(valid_idx) = mafdr(p_WW(valid_idx), 'BHFDR', true);
pFDR_WW = reshape(p_tmp, size(p_WW));

%% ---------- Part 2: G-W permutation（100×68 → 7×5） ----------
fprintf('Permutation for G-W blocks...\n');

Net_perm_GW = nan(7, 5, nPerm);

for i_perm = 1:nPerm
    Net_sets = nan(7, 5, nSet);
    
    for i_set = 1:nSet
        current_set = datasets_str{i_set};
        data_dir    = datasets_dir{i_set};
        this_perm_dir = perm_dirs_all{i_set}{i_perm};

        fold_files = g_ls([this_perm_dir '/GWFC/Fold_*_Score.mat']);
        clear W_cv
        for i_cv = 1:length(fold_files)
            tmp = load(fold_files{i_cv});
            W_cv(i_cv, :) = tmp.w_Brain_Haufe(:)'; 
        end
        w_perm = mean(W_cv, 1);

        % vector -> 100x68
        nG = 100; nW = 68;
        if numel(w_perm) ~= nG*nW
            error('GW vector length mismatch in %s perm %d', current_set, i_perm);
        end

        if strcmpi(current_set, 'EFNY')
           
            FC_gw = reshape(w_perm, [nW, nG])';
        else
            
            FC_gw = reshape(w_perm, [nG, nW]);
        end

        % 100x68 → 7x5 
        Net_tmp = zeros(7,5);
        for g_row = 1:7
            for w_col = 1:5
                idx_g = (gm_labels == g_row);   % 100
                idx_w = (wm_labels == w_col);   % 68
                subM  = FC_gw(idx_g, idx_w);
                Net_tmp(g_row, w_col) = mean(subM(:));
            end
        end

        Net_sets(:,:,i_set) = Net_tmp;
    end

    Net_perm_GW(:,:,i_perm) = mean(Net_sets, 3);
end

% cal G-W p value
p_GW    = nan(7,5);
pFDR_GW = nan(7,5);
for i = 1:7
    for j = 1:5
        real_val  = Net_real_GW(i,j);
        null_vals = squeeze(Net_perm_GW(i,j,:));

        if real_val > 0
            p_GW(i,j) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_GW(i,j) = sum(null_vals <= real_val) / nPerm;
        else
            p_GW(i,j) = NaN;
        end
    end
end
valid_idx = ~isnan(p_GW(:));
p_tmp = nan(size(p_GW));
p_tmp(valid_idx) = mafdr(p_GW(valid_idx), 'BHFDR', true);
pFDR_GW = reshape(p_tmp, size(p_GW));

%% ---------- Part 3: G-G permutation（100×100 → 7×7） ----------
fprintf('Permutation for G-G blocks...\n');

Net_perm_GG = nan(7, 7, nPerm);

for i_perm = 1:nPerm
    Net_sets = nan(7, 7, nSet);
    
    for i_set = 1:nSet
        current_set = datasets_str{i_set};
        data_dir    = datasets_dir{i_set};
        this_perm_dir = perm_dirs_all{i_set}{i_perm};

        fold_files = g_ls([this_perm_dir '/GGFC/Fold_*_Score.mat']);
        clear W_cv
        for i_cv = 1:length(fold_files)
            tmp = load(fold_files{i_cv});
            W_cv(i_cv, :) = tmp.w_Brain_Haufe(:)';  
        end
        w_perm = mean(W_cv, 1);

        % vector -> 100x100
        m = 100; n = 100;
        if numel(w_perm) ~= m*(m-1)/2
            error('GG vector length mismatch in %s perm %d', current_set, i_perm);
        end

        FC = zeros(m, n);
        if strcmpi(current_set, 'EFNY')
           
            [ri, ci] = find(tril(true(m, n), -1));
            [~, ord] = sortrows([ri, ci], [1 2]);
            lin_idx   = sub2ind([m n], ri(ord), ci(ord));
            FC(lin_idx) = w_perm(:);
        else
           
            tril_idx = find(tril(true(m, n), -1));
            FC(tril_idx) = w_perm(:);
        end
        FC_full = FC + FC';
        FC_full(logical(eye(m))) = 0;

        % 100x100 → 7x7 
        Net_tmp = zeros(7,7);
        for g_row = 1:7
            for g_col = 1:7
                idx_r = (gm_labels == g_row);
                idx_c = (gm_labels == g_col);
                subM  = FC_full(idx_r, idx_c);
                Net_tmp(g_row, g_col) = mean(subM(:));
            end
        end

        Net_sets(:,:,i_set) = Net_tmp;
    end

    Net_perm_GG(:,:,i_perm) = mean(Net_sets, 3);
end

% cal G-G 的 permutation p value
p_GG    = nan(7,7);
pFDR_GG = nan(7,7);
for i = 1:7
    for j = 1:7
        real_val  = Net_real_GG(i,j);
        null_vals = squeeze(Net_perm_GG(i,j,:));

        if real_val > 0
            p_GG(i,j) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_GG(i,j) = sum(null_vals <= real_val) / nPerm;
        else
            p_GG(i,j) = NaN;
        end
    end
end
valid_idx = ~isnan(p_GG(:));
p_tmp = nan(size(p_GG));
p_tmp(valid_idx) = mafdr(p_GG(valid_idx), 'BHFDR', true);
pFDR_GG = reshape(p_tmp, size(p_GG));

%% ---------- save ----------
save([out_dir '/datasetsAveraged_perm_blockStats_WW_GW_GG_mean.mat'], ...
    'Net_real_WW', 'Net_perm_WW', 'p_WW', 'pFDR_WW', ...
    'Net_real_GW', 'Net_perm_GW', 'p_GW', 'pFDR_GW', ...
    'Net_real_GG', 'Net_perm_GG', 'p_GG', 'pFDR_GG', ...
    'nPerm', 'datasets_str');