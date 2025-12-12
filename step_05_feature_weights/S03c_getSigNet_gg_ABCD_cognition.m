%% ABCD cognition: G-G network-level permutation test
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
addpath([proj_dir '/ABCD/code/6th_post_prediction']);

current_set = 'ABCDcog';
cog_str     = 'nihtbx_totalcomp';
data_dir    = [proj_dir '/ABCD/code/5th_prediction/PLS_prediction/networkFC/' ...
               'cognition_nosmooth_motion2runFD/' cog_str];
weight_dir  = [proj_dir '/ABCD/results/figure/weight_cognition'];
real_mat_file = [weight_dir '/' current_set '_Haufe_FCmatrix_gg_Schaefer100.mat'];

out_dir     = [proj_dir '/ABCD/results/figure/weight_cognition'];
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% ---------- 0. 基本信息 & 标签 ----------
m = 100; n = 100;
nEdge = m * (m - 1) / 2;

load gm_labels.mat;           % 100x1, 1~7
load cmap0.mat;
load Schaefer100_info.mat;    % networkName: 7 GM 
if numel(gm_labels) ~= m
    error('gm_labels length (%d) != 100', numel(gm_labels));
end
nNet = 7;

%% ---------- 1. real G-G FC matrix ----------
VisualizeFolder = [data_dir '/WeightVisualize_Overall_RandomCV'];
load([VisualizeFolder '/w_Brain_median_101CV_GG_Haufe.mat']);   % w_Brain_median_101CV

w_obs_vec = w_Brain_median_101CV(:);
if numel(w_obs_vec) ~= nEdge
    error('GG The number of edges is incorrect: nEdge = %d, but 100×99/2 = %d', numel(w_obs_vec), nEdge);
end


[ri, ci] = find(tril(true(m, n), -1));     
[~, ord] = sortrows([ri, ci], [1 2]);   
ri = ri(ord); ci = ci(ord);
lin_idx = sub2ind([m, n], ri, ci);

FC = zeros(m, n);
FC(lin_idx) = w_obs_vec(:);    
FC_real = FC + FC';
FC_real(logical(eye(m))) = 0;   % 100x100 

Net_real_GG = zeros(nNet, nNet);
for g_row = 1:nNet
    idx_r = (gm_labels == g_row);
    for g_col = 1:nNet
        idx_c = (gm_labels == g_col);
        subM  = FC_real(idx_r, idx_c);
        % Net_real_GG(g_row, g_col) = mean(subM(:), 'omitnan');
        Net_real_GG(g_row, g_col) = mean(subM(:));
    end
end

%% ---------- 2. permutation ----------
Time_dirs = g_ls([data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
nPerm     = length(Time_dirs);
if nPerm == 0
    error('Time_* directory not found: %s', ...
          [data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
end
fprintf('GG: %d permutations (Time_*) detected.\n', nPerm);

Net_perm_GG = nan(nNet, nNet, nPerm);

for i_perm = 1:nPerm
    this_time_dir = Time_dirs{i_perm};
    fold_files = g_ls([this_time_dir '/GGFC/Fold_*_Score.mat']);
    if isempty(fold_files)
        error('The file GGFC/Fold_*_Score.mat is not present under Time %d (%s).', ...
              i_perm, this_time_dir);
    end

    clear W_cv
    for i_fold = 1:length(fold_files)
        tmp = load(fold_files{i_fold});
        if ~isfield(tmp, 'w_Brain_Haufe')
            error('The file %s does not contain w_Brain_Haufe', fold_files{i_fold});
        end
        this_w = tmp.w_Brain_Haufe(:)';
        if numel(this_w) ~= nEdge
            error('Time %d, fold %d: edge count mismatch (%d vs %d)', ...
                  i_perm, i_fold, numel(this_w), nEdge);
        end
        W_cv(i_fold, :) = this_w;  
    end

    w_perm = mean(W_cv, 1);        % 1×nEdge

    % vector -> 100x100
    FC = zeros(m, n);
    FC(lin_idx) = w_perm(:);
    FC_full = FC + FC';
    FC_full(logical(eye(m))) = 0;

    % 100x100 → 7x7
    Net_tmp = zeros(nNet, nNet);
    for g_row = 1:nNet
        idx_r = (gm_labels == g_row);
        for g_col = 1:nNet
            idx_c = (gm_labels == g_col);
            subM  = FC_full(idx_r, idx_c);
            % Net_tmp(g_row, g_col) = mean(subM(:), 'omitnan');
            Net_tmp(g_row, g_col) = mean(subM(:));
        end
    end

    Net_perm_GG(:, :, i_perm) = Net_tmp;
end

%% ---------- 3. calculate G-G  permutation p value ----------
p_GG    = nan(nNet, nNet);
pFDR_GG = nan(nNet, nNet);

for i = 1:nNet
    for j = 1:nNet
        real_val  = Net_real_GG(i, j);
        null_vals = squeeze(Net_perm_GG(i, j, :));

        if all(isnan(null_vals))
            p_GG(i, j) = NaN;
            continue;
        end

        if real_val > 0
            p_GG(i, j) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_GG(i, j) = sum(null_vals <= real_val) / nPerm;
        else
            p_GG(i, j) = NaN;
        end
    end
end

valid_mask = ~isnan(p_GG);
p_vec      = p_GG(valid_mask);
p_vec_fdr  = mafdr(p_vec, 'BHFDR', true);
pFDR_GG(valid_mask) = p_vec_fdr;

alpha_raw   = 0.05;
Net_sig_raw = (p_GG <= alpha_raw) & valid_mask;

[row_sig, col_sig] = find(tril(Net_sig_raw));
fprintf('==== GG Network-level significant pairs (raw p < %.2f) ====\n', alpha_raw);
for k = 1:numel(row_sig)
    i = row_sig(k);
    j = col_sig(k);
    fprintf('GM Net %d (%s) - GM Net %d (%s): mean = %.4f, p = %.4g, p_FDR = %.4g\n', ...
        i, networkName{i}, j, networkName{j}, ...
        Net_real_GG(i,j), p_GG(i,j), pFDR_GG(i,j));
end

% save
save_name_net = [out_dir '/' current_set '_Haufe_GG_networkLevel_perm_Schaefer100_ageStyle'];
save(save_name_net, 'Net_real_GG', 'Net_perm_GG', ...
    'p_GG', 'pFDR_GG', 'Net_sig_raw', 'nPerm');

%% ---------- 4.plot 7x7 GM ----------
figure('Position', [150, 150, 450, 400]);
imagesc(Net_real_GG);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';

maxAbs = max(abs(Net_real_GG(~isnan(Net_real_GG))));
if isempty(maxAbs) || maxAbs == 0
    maxAbs = 1e-4;
end
caxis([-maxAbs, maxAbs]);
colormap(cmap);

xticks(1:nNet);
xticklabels(networkName);
yticks(1:nNet);
yticklabels(networkName);

ax = gca;
ax.XTickLabelRotation = 45;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title(['G-G network-level Haufe weights of ' current_set]);
box off;
axis square;

%*: raw p < 0.05
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

saveas(gcf, [save_name_net '_7x7.fig']);
saveas(gcf, [save_name_net '_7x7.jpg']);
print(gcf,  [save_name_net '_7x7.pdf'], '-dpdf', '-bestfit');