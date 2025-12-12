%% ABCD pfactor: G-G network-level permutation test
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
addpath([proj_dir '/ABCD/code/6th_post_prediction']);

current_set = 'ABCDpfactor';
cog_str = 'ADHD';
data_dir = [proj_dir '/ABCD/code/5th_prediction/PLS_prediction/networkFC/p_factor_nosmooth_motion2runFD/' cog_str];

weight_dir  = [proj_dir '/ABCD/results/figure/weight_pfactor'];
real_mat_file = [weight_dir '/' current_set '_Haufe_FCmatrix_ww_Schaefer100.mat'];

out_dir = [proj_dir '/ABCD/results/figure/weight_pfactor'];
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ============ 0. label ============
load wm_labels.mat;           % 68x1, range 1~5
load cmap0.mat;
load JHU68_info_rename.mat;   % tractsName
nNet = 5;

m = 68; n = 68;

%% ============ 1. real 5x5 matrix ============
S = load(real_mat_file, 'FCmatrix_full');
FC_real = S.FCmatrix_full;   % 68x68

if ~isequal(size(FC_real), [m, n])
    error('The size of FCmatrix_full is not 68x68, but %dx%d', size(FC_real,1), size(FC_real,2));
end

Net_real_WW = zeros(nNet, nNet);
for net1 = 1:nNet
    idx_r = (wm_labels == net1);
    for net2 = 1:nNet
        idx_c = (wm_labels == net2);
        subM  = FC_real(idx_r, idx_c);
        Net_real_WW(net1, net2) = mean(subM(:));
    end
end

%% ============ 2. vector --> matrix ============
[ri, ci] = find(tril(true(m, n), -1));      
[~, ord] = sortrows([ri, ci], [1 2]);       
ri = ri(ord);
ci = ci(ord);
lin_idx = sub2ind([m, n], ri, ci);
nEdge = numel(lin_idx);

%% ============ 3. read permutation（Time_*） ============
Time_dirs = g_ls([data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
nPerm     = length(Time_dirs);
if nPerm == 0
    error('The Time_* directory was not found. Please check the path: %s', ...
          [data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
end
fprintf('%d permutations (Time_*) have been detected', nPerm);

Net_perm_WW = nan(nNet, nNet, nPerm);

for i_perm = 1:nPerm
    this_time_dir = Time_dirs{i_perm};

    fold_files = g_ls([this_time_dir '/WWFC/Fold_*_Score.mat']);
    if isempty(fold_files)
        error('The file WWFC/Fold_*_Score.mat was not found under Time %d (%s)', ...
              i_perm, this_time_dir);
    end

    clear W_cv
    for i_fold = 1:length(fold_files)
        tmp = load(fold_files{i_fold});
        if ~isfield(tmp, 'w_Brain_Haufe')
            error('The variable w_Brain_Haufe is not find in the %s file', fold_files{i_fold});
        end
        this_w = tmp.w_Brain_Haufe(:)';
        if i_fold == 1
            if numel(this_w) ~= nEdge
                error('Time %d, fold %d: The number of edges is incorrect (%d vs %d).', ...
                      i_perm, i_fold, numel(this_w), nEdge);
            end
        else
            if numel(this_w) ~= nEdge
                error('Time %d, fold %d: The number of sides is inconsistent (%d vs %d)', ...
                      i_perm, i_fold, numel(this_w), nEdge);
            end
        end
        W_cv(i_fold, :) = this_w;   
    end

    w_perm = mean(W_cv, 1);   % 1 × nEdge

    % vector -> 68x68
    FC = zeros(m, n);
    FC(lin_idx) = w_perm(:);
    FC_full = FC + FC';
    FC_full(logical(eye(m))) = 0;

    % 68x68 → 5x5 
    Net_tmp = zeros(nNet, nNet);
    for net1 = 1:nNet
        idx_r = (wm_labels == net1);
        for net2 = 1:nNet
            idx_c = (wm_labels == net2);
            subM  = FC_full(idx_r, idx_c);
            Net_tmp(net1, net2) = mean(subM(:));
        end
    end

    Net_perm_WW(:, :, i_perm) = Net_tmp;
end

%% ============ 4. calculate permutation p value ============
p_WW    = nan(nNet, nNet);
pFDR_WW = nan(nNet, nNet);

for i = 1:nNet
    for j = 1:nNet
        real_val  = Net_real_WW(i, j);
        null_vals = squeeze(Net_perm_WW(i, j, :));

        if all(isnan(null_vals))
            p_WW(i, j) = NaN;
            continue;
        end

        if real_val > 0
            p_WW(i, j) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_WW(i, j) = sum(null_vals <= real_val) / nPerm;
        else
            p_WW(i, j) = NaN;
        end
    end
end

% FDR corrected
valid_mask = ~isnan(p_WW);
p_vec      = p_WW(valid_mask);
p_vec_fdr  = mafdr(p_vec, 'BHFDR', true);
pFDR_WW(valid_mask) = p_vec_fdr;

alpha_raw = 0.05;
Net_sig_raw = (p_WW <= alpha_raw) & valid_mask;   % raw p < 0.05

% print significance
[row_sig, col_sig] = find(tril(Net_sig_raw));
fprintf('==== Network-level significant pairs (raw p < %.2f) ====\n', alpha_raw);
for k = 1:numel(row_sig)
    i = row_sig(k);
    j = col_sig(k);
    fprintf('Net %d (%s) - Net %d (%s): mean = %.4f, p = %.4g, p_FDR = %.4g\n', ...
            i, tractsName{i}, j, tractsName{j}, ...
            Net_real_WW(i,j), p_WW(i,j), pFDR_WW(i,j));
end

% save
save_name_net = [out_dir '/' current_set '_Haufe_WW_networkLevel_perm_Schaefer100_ageStyle'];
save(save_name_net, 'Net_real_WW', 'Net_perm_WW', 'p_WW', 'pFDR_WW', ...
                    'Net_sig_raw', 'nPerm');

%% ============ 5. plot 5×5 ============
figure('Position', [150, 150, 450, 400]);
imagesc(Net_real_WW);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';

maxAbs = max(abs(Net_real_WW(~isnan(Net_real_WW))));
if isempty(maxAbs) || maxAbs == 0
    maxAbs = 1e-4;
end
caxis([-maxAbs, maxAbs]);
colormap(cmap);

xticks(1:nNet);
xticklabels(tractsName);
yticks(1:nNet);
yticklabels(tractsName);

ax = gca;
ax.XTickLabelRotation = 0;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title(['W-W network-level Haufe weights of ' current_set]);
box off;
axis square;

% *：raw p < 0.05
hold on;
for i = 1:nNet
    for j = 1:nNet
        if Net_sig_raw(i, j)
            text(j, i, '*', 'HorizontalAlignment', 'center', ...
                       'VerticalAlignment', 'middle', ...
                       'FontSize', 18, 'FontWeight', 'bold', ...
                       'Color', 'k');
        end
    end
end
hold off;

saveas(gcf, [save_name_net '_5x5.fig']);
saveas(gcf, [save_name_net '_5x5.jpg']);
print(gcf,  [save_name_net '_5x5.pdf'], '-dpdf', '-bestfit');