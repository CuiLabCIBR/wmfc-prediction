%% ===============================================================
%  ABCD pfactor (ADHD):W-W / G-W / G-G
%% ===============================================================
clear;
proj_dir    = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
addpath([proj_dir '/ABCD/code/6th_post_prediction']);  % Schaefer100_info / JHU68_info_rename 等

current_set = 'ABCDpfactor';

% permutation results dir
res_dir = [proj_dir '/ABCD/results/figure/weight_pfactor'];

% output path
fig_dir = [res_dir '/NetBlock_plots'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% ---------- 1. load pfactor–ADHD permutation results ----------
file_WW = [res_dir '/' current_set '_Haufe_WW_networkLevel_perm_Schaefer100_ageStyle.mat'];
file_GW = [res_dir '/' current_set '_Haufe_GW_networkLevel_perm_Schaefer100_ageStyle.mat'];
file_GG = [res_dir '/' current_set '_Haufe_GG_networkLevel_perm_Schaefer100_ageStyle.mat'];

S_WW = load(file_WW, 'Net_real_WW', 'p_WW', 'pFDR_WW');
S_GW = load(file_GW, 'Net_real_GW', 'p_GW', 'pFDR_GW');
S_GG = load(file_GG, 'Net_real_GG', 'p_GG', 'pFDR_GG');

%% ---------- 2. load colormap  ----------
load cmap0.mat;   

%% ---------- 3. load GM / WM netname ----------
load Schaefer100_info.mat;    % networkName: 7 GM
load JHU68_info_rename.mat;   % tractsName: 5 WM 

GM_net_names = networkName;  
WM_net_names = tractsName;   

%% ---------- 4. plot(pfactor–ADHD) ----------
% ---- 4.1 W-W: 5×5 ----
saveName_WW = [fig_dir '/ABCDpfactor_ADHD_NetBlock_WW_weight'];
plot_block_matrix(S_WW.Net_real_WW, S_WW.p_WW, ...
                  WM_net_names, WM_net_names, ...
                  'full', 'ABCD pfactor (ADHD) W-W', saveName_WW, cmap);

% ---- 4.2 G-W: 7×5 ----
saveName_GW = [fig_dir '/ABCDpfactor_ADHD_NetBlock_GW_weight'];
plot_block_matrix(S_GW.Net_real_GW, S_GW.p_GW, ...
                  GM_net_names, WM_net_names, ...
                  'full', 'ABCD pfactor (ADHD) G-W', saveName_GW, cmap);

% ---- 4.3 G-G: 7×7 ----
saveName_GG = [fig_dir '/ABCDpfactor_ADHD_NetBlock_GG_weight'];
plot_block_matrix(S_GG.Net_real_GG, S_GG.p_GG, ...
                  GM_net_names, GM_net_names, ...
                  'full', 'ABCD pfactor (ADHD) G-G', saveName_GG, cmap);