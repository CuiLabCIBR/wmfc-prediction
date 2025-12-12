%% ===============================================================
%  ABCD cognition: W-W / G-W / G-G
%% ===============================================================
clear;
proj_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath([proj_dir '/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64']);
addpath([proj_dir '/ABCD/code/6th_post_prediction']);  

current_set = 'ABCDcog';


res_dir = [proj_dir '/ABCD/results/figure/weight_cognition'];


fig_dir = [res_dir '/NetBlock_plots'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% ---------- 1. load ABCD permutation results ----------
file_WW = [res_dir '/' current_set '_Haufe_WW_networkLevel_perm_Schaefer100_ageStyle.mat'];
file_GW = [res_dir '/' current_set '_Haufe_GW_networkLevel_perm_Schaefer100_ageStyle.mat'];
file_GG = [res_dir '/' current_set '_Haufe_GG_networkLevel_perm_Schaefer100_ageStyle.mat'];

S_WW = load(file_WW, 'Net_real_WW', 'p_WW');
S_GW = load(file_GW, 'Net_real_GW', 'p_GW');
S_GG = load(file_GG, 'Net_real_GG', 'p_GG');

%% ---------- 2. load colormap ----------
load cmap0.mat;   

%% ---------- 3. load GM / WM name ----------
load Schaefer100_info.mat;
load JHU68_info_rename.mat;  

GM_net_names = networkName; 
WM_net_names = tractsName;   

%% ---------- 4. plot ---------
% ---- W-W: 5×5， ----
saveName_WW = [fig_dir '/ABCDcog_NetBlock_WW_weight'];
plot_block_matrix(S_WW.Net_real_WW, S_WW.p_WW, ...
                  WM_net_names, WM_net_names, ...
                  'full', 'ABCD W-W', saveName_WW, cmap);

% ---- G-W: 7×5----
saveName_GW = [fig_dir '/ABCDcog_NetBlock_GW_weight'];
plot_block_matrix(S_GW.Net_real_GW, S_GW.p_GW, ...
                  GM_net_names, WM_net_names, ...
                  'full', 'ABCD G-W', saveName_GW, cmap);

% ---- G-G: 7×7 ----
saveName_GG = [fig_dir '/ABCDcog_NetBlock_GG_weight'];
plot_block_matrix(S_GG.Net_real_GG, S_GG.p_GG, ...
                  GM_net_names, GM_net_names, ...
                  'full', 'ABCD G-G', saveName_GG, cmap);