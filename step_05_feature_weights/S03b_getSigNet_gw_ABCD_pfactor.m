%% ABCD pfactor: G-W network-level permutation test 
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
addpath([proj_dir '/ABCD/code/6th_post_prediction']);

current_set = 'ABCDpfactor';
cog_str = 'ADHD';
data_dir = [proj_dir '/ABCD/code/5th_prediction/PLS_prediction/networkFC/p_factor_nosmooth_motion2runFD/' cog_str];

weight_dir  = [proj_dir '/ABCD/results/figure/weight_pfactor'];
real_mat_file = [weight_dir '/' current_set '_Haufe_FCmatrix_gw_Schaefer100.mat'];

out_dir = [proj_dir '/ABCD/results/figure/weight_pfactor'];
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end


%% ---------- 0. label ----------
m_gm = 100;          % Schaefer100
n_wm = 68;           % JHU68
nEdge = m_gm * n_wm;

load gm_labels.mat;           % 100x1, 1~7
load wm_labels.mat;           % 68x1, 1~5
load cmap0.mat;
load Schaefer100_info.mat;    % networkName
load JHU68_info_rename.mat;   % tractsName

if numel(gm_labels) ~= m_gm
    error('gm_labels length(%d) != 100', numel(gm_labels));
end
if numel(wm_labels) ~= n_wm
    error('wm_labels length(%d) != 68', numel(wm_labels));
end

nGM = 7;
nWM = 5;

%% ---------- 1. real G-W matrix: Net_real_GW ----------
VisualizeFolder = [data_dir '/WeightVisualize_Overall_RandomCV'];
load([VisualizeFolder '/w_Brain_median_101CV_GW_Haufe.mat']);   % w_Brain_median_101CV

w_obs_vec = w_Brain_median_101CV(:);
if numel(w_obs_vec) ~= nEdge
    error('The edge count in GW is incorrect: nEdge = %d, but 100×68 = %d', numel(w_obs_vec), nEdge);
end

% Python row-major: vec = FC(1,1..68), FC(2,1..68), ...
FC_real = reshape(w_obs_vec, [n_wm, m_gm])';   % 100x68

Net_real_GW = zeros(nGM, nWM);
for g_net = 1:nGM
    idx_g = (gm_labels == g_net);
    for w_net = 1:nWM
        idx_w = (wm_labels == w_net);
        subM  = FC_real(idx_g, idx_w);
        Net_real_GW(g_net, w_net) = mean(subM(:));
    end
end

%% ---------- 2. permutation ----------
Time_dirs = g_ls([data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
nPerm     = length(Time_dirs);
if nPerm == 0
    error('The Time_* directory was found: %s', ...
          [data_dir '/RegressCovariates_RandomCV_Permutation/Time_*']);
end
fprintf('GW: find %d  permutation (Time_*).\n', nPerm);

Net_perm_GW = nan(nGM, nWM, nPerm);   % 7x5xNperm

for i_perm = 1:nPerm
    this_time_dir = Time_dirs{i_perm};
    fold_files = g_ls([this_time_dir '/GWFC/Fold_*_Score.mat']);
    if isempty(fold_files)
        error('There is no GWFC/Fold_*_Score.mat file under Time %d (%s).', ...
              i_perm, this_time_dir);
    end

    clear W_cv
    for i_fold = 1:length(fold_files)
        tmp    = load(fold_files{i_fold});
        if ~isfield(tmp, 'w_Brain_Haufe')
            error('The file %s does not contain w_Brain_Haufe', fold_files{i_fold});
        end
        this_w = tmp.w_Brain_Haufe(:)';
        if numel(this_w) ~= nEdge
            error('Time %d, fold %d: The number of edges is inconsistent (%d vs %d)', ...
                  i_perm, i_fold, numel(this_w), nEdge);
        end
        W_cv(i_fold, :) = this_w;    
    end


    w_perm = mean(W_cv, 1);         % 1×nEdge

    % vector -> 100x68
    FC_perm = reshape(w_perm, [n_wm, m_gm])';   % 100x68

    % 100x68 → 7x5 
    Net_tmp = zeros(nGM, nWM);
    for g_net = 1:nGM
        idx_g = (gm_labels == g_net);
        for w_net = 1:nWM
            idx_w = (wm_labels == w_net);
            subM  = FC_perm(idx_g, idx_w);
            Net_tmp(g_net, w_net) = mean(subM(:));
        end
    end

    Net_perm_GW(:, :, i_perm) = Net_tmp;
end

%% ---------- 3. calculate permutation p value ----------
p_GW    = nan(nGM, nWM);
pFDR_GW = nan(nGM, nWM);

for gi = 1:nGM
    for wi = 1:nWM
        real_val  = Net_real_GW(gi, wi);
        null_vals = squeeze(Net_perm_GW(gi, wi, :));

        if all(isnan(null_vals))
            p_GW(gi, wi) = NaN;
            continue;
        end

        if real_val > 0
            p_GW(gi, wi) = sum(null_vals >= real_val) / nPerm;
        elseif real_val < 0
            p_GW(gi, wi) = sum(null_vals <= real_val) / nPerm;
        else
            p_GW(gi, wi) = NaN;
        end
    end
end

valid_mask = ~isnan(p_GW);
p_vec      = p_GW(valid_mask);
p_vec_fdr  = mafdr(p_vec, 'BHFDR', true);
pFDR_GW(valid_mask) = p_vec_fdr;

alpha_raw   = 0.05;
Net_sig_raw = (p_GW <= alpha_raw) & valid_mask;

% print significance
[g_sig, w_sig] = find(Net_sig_raw);
fprintf('==== GW Network-level significant pairs (raw p < %.2f) ====\n', alpha_raw);
for k = 1:numel(g_sig)
    gi = g_sig(k);
    wi = w_sig(k);
    fprintf('GM Net %d (%s) - WM Net %d (%s): mean = %.4f, p = %.4g, p_FDR = %.4g\n', ...
        gi, networkName{gi}, wi, tractsName{wi}, ...
        Net_real_GW(gi,wi), p_GW(gi,wi), pFDR_GW(gi,wi));
end

% save
save_name_net = [out_dir '/' current_set '_Haufe_GW_networkLevel_perm_Schaefer100_ageStyle'];
save(save_name_net, 'Net_real_GW', 'Net_perm_GW', ...
    'p_GW', 'pFDR_GW', 'Net_sig_raw', 'nPerm');

%% ---------- 4. plot 7x5 GM×WM  ----------
figure('Position', [200, 150, 520, 420]);
imagesc(Net_real_GW);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';

maxAbs = max(abs(Net_real_GW(~isnan(Net_real_GW))));
if isempty(maxAbs) || maxAbs == 0
    maxAbs = 1e-4;
end
caxis([-maxAbs, maxAbs]);
colormap(cmap);

xticks(1:nWM);
xticklabels(tractsName);
yticks(1:nGM);
yticklabels(networkName);

ax = gca;
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title(['G-W network-level Haufe weights of ' current_set]);
box off;
axis tight; axis ij;

% *: raw p < 0.05
hold on;
for gi = 1:nGM
    for wi = 1:nWM
        if Net_sig_raw(gi, wi)
            text(wi, gi, '*', 'HorizontalAlignment', 'center', ...
                          'VerticalAlignment', 'middle', ...
                          'FontSize', 18, 'FontWeight', 'bold', ...
                          'Color', 'k');
        end
    end
end
hold off;

saveas(gcf, [save_name_net '_GMxWM_7x5.fig']);
saveas(gcf, [save_name_net '_GMxWM_7x5.jpg']);
print(gcf,  [save_name_net '_GMxWM_7x5.pdf'], '-dpdf', '-bestfit');