%% ===============================================================
%  plot age feature weights（W-W / G-W / G-G）
%% ===============================================================
clear;
proj_dir = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');

out_dir  = [proj_dir '/Figure/Age_perdiction_feature_weight'];

%% ---------- 1. load permutation results ----------
permFile = [out_dir '/datasetsAveraged_perm_blockStats_WW_GW_GG_mean.mat'];
load(permFile, ...
    'Net_real_WW', 'Net_real_GW', 'Net_real_GG', ...
    'p_WW', 'p_GW', 'p_GG');

%% ---------- 2.load colormap ----------

load cmap0.mat;

%% ---------- 3. load label and network name ----------

WM_net_names = {'BS','COMM','ASSOC','PROJ','SWM'};
GM_net_names = {'VIS','SMN','DAN','VAN','LIM','CON','DMN'};

%% ---------- 5. plot ----------
% W-W：5×5
saveName_WW = [out_dir '/NetBlock_WW_weight'];
plot_block_matrix(Net_real_WW, p_WW, WM_net_names, WM_net_names, ...
                  'full', 'W-W', saveName_WW, cmap);

% G-W：7×5
saveName_GW = [out_dir '/NetBlock_GW_weight'];
plot_block_matrix(Net_real_GW, p_GW, GM_net_names, WM_net_names, ...
                  'full', 'G-W', saveName_GW, cmap);

% G-G：7×7
saveName_GG = [out_dir '/NetBlock_GG_weight'];
% plot_block_matrix(Net_real_GG, p_GG, GM_net_names, GM_net_names, ...
%                   'lower', 'G-G', saveName_GG, cmap);
plot_block_matrix(Net_real_GG, p_GG, GM_net_names, GM_net_names, ...
                  'full', 'G-G', saveName_GG, cmap);

