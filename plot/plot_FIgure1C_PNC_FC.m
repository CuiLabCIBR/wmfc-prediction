%% PNC FC matrices: GG / GW / WW (group mean)
clear; clc; close all;

load subid_final_PNC.mat;
subList = subid_final_PNC;

GM_network_num = 100; 
WM_network_num = 68; 

total_FCmatrix_WW = zeros(WM_network_num,WM_network_num,length(subList));
total_FCmatrix_GW = zeros(GM_network_num,WM_network_num,length(subList));
total_FCmatrix_GG = zeros(GM_network_num,GM_network_num,length(subList));

datapath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/results';
savepath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/Figure/PNC_FC';
%%
for i_sub=1:length(subList)
    i_sub
    sub_current = subList{i_sub};
    %% 1. WM-WM
    %load([datapath '/FCmatrix_individual/' sub_current '/' sub_current '_rICBM_region50_WW_FC.mat']);
	load([datapath '/FCmatrix_individual_schaefer100_nosmooth/' sub_current '/' sub_current '_rICBM_60p_WW_FC.mat']);
    WW_FC = WW_FCmatrix;
    WW_FC(isnan(WW_FC)) = 0;
    total_FCmatrix_WW(:,:,i_sub) = WW_FC;
       
    %% 2. GM-WM
    %load([datapath '/FCmatrix_individual/' sub_current '/' sub_current '_Yeo7_rICBM_GW_FC.mat']);
	load([datapath '/FCmatrix_individual_schaefer100_nosmooth/' sub_current '/' sub_current '_Schaefer100_rICBM_60p_GW_FC.mat']);
    GW_FC = GW_FCmatrix;
    GW_FC(isnan(GW_FC)) = 0;
    total_FCmatrix_GW(:,:,i_sub) = GW_FC;
    
   %% 3. GM-GM
    load([datapath '/FCmatrix_individual_schaefer100_nosmooth/' sub_current '/' sub_current '_Schaefer100_GG_FC.mat']);
    GG_FC = GG_FCmatrix;
    GG_FC(isnan(GG_FC)) = 0;
    total_FCmatrix_GG(:,:,i_sub) = GG_FC;
    
end


WW = mean(total_FCmatrix_WW, 3);


GW = mean(total_FCmatrix_GW, 3);

GG = mean(total_FCmatrix_GG, 3);


GG(1:size(GG,1)+1:end) = 0;
WW(1:size(WW,1)+1:end) = 0;


load cmap0.mat                   
load Schaefer100_info.mat       
load JHU68_info_rename.mat     

%% ======================================================================
%                               1. G-G
% =======================================================================
FC_GG = GG(regionID_sortedByNetwork, regionID_sortedByNetwork);

figure('Position', [100, 100, 550, 500]);
imagesc(FC_GG);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';
colormap(cmap);
caxis([-max(abs(FC_GG(:))), max(abs(FC_GG(:)))]);

hold on;
% Yeo7 
split_positions_GG = [17, 31, 46, 58, 63, 76];
for i = 1:length(split_positions_GG)
    pos = split_positions_GG(i) + 0.5;

    line([0.5, 100.5], [pos, pos], 'Color', 'k', 'LineWidth', 0.5);

    line([pos, pos], [0.5, 100.5], 'Color', 'k', 'LineWidth', 0.5);
end
hold off;

% network label
xtick_positions_GG = [8, 25, 39, 52, 61, 70, 90];
ytick_positions_GG = [8, 25, 39, 52, 61, 70, 90];

xticks(xtick_positions_GG);
yticks(ytick_positions_GG);
xticklabels(networkName);
yticklabels(networkName);

ax = gca;
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title('G-G FC');

save_name_GG = fullfile(savepath, 'PNC_GG_FCmean');
saveas(gcf, [save_name_GG '.fig']);
saveas(gcf, [save_name_GG '.jpg']);
print(gcf, [save_name_GG '.pdf'], '-dpdf', '-bestfit');

%% ======================================================================
%                               2. G-W
% =======================================================================
%GW : 100x68（row=GM，column=WM）：
FC_GW = GW(regionID_sortedByNetwork, regionID_sortedByTracts);

vals_GW   = FC_GW(:);
num_pos   = sum(vals_GW > 0);
num_neg   = sum(vals_GW < 0);
num_zero  = sum(vals_GW == 0);
num_total = numel(vals_GW);

fprintf('G-W FC ：\n');
fprintf('  positive: %d (%.2f%%)\n',  num_pos,  100 * num_pos  / num_total);
fprintf('  negative: %d (%.2f%%)\n',  num_neg,  100 * num_neg  / num_total);
fprintf('  zero: %d (%.2f%%)\n\n', num_zero, 100 * num_zero / num_total);

figure('Position', [100, 100, 550, 500]);
imagesc(FC_GW);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';
colormap(cmap);
caxis([-max(abs(FC_GW(:))), max(abs(FC_GW(:)))]);

hold on;

col_separators = [10, 15, 34, 50];  
for i = 1:length(col_separators)
    x = col_separators(i) + 0.5;
    line([x, x], [0.5, 100.5], 'Color', 'k', 'LineWidth', 0.5);
end


row_separators = [17, 31, 46, 58, 63, 76];
for i = 1:length(row_separators)
    y = row_separators(i) + 0.5;
    line([0.5, 68.5], [y, y], 'Color', 'k', 'LineWidth', 0.5);
end
hold off;


xtick_positions_GW = [5, 13, 25, 43, 59];
xticks(xtick_positions_GW);
xticklabels(tractsName);


ytick_positions_GW = [8, 25, 39, 52, 61, 70, 90];
yticks(ytick_positions_GW);
yticklabels(networkName);

ax = gca;
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title('G-W FC');

save_name_GW = fullfile(savepath, 'PNC_GW_FCmean');
saveas(gcf, [save_name_GW '.fig']);
saveas(gcf, [save_name_GW '.jpg']);
print(gcf, [save_name_GW '.pdf'], '-dpdf', '-bestfit');

%% ======================================================================
%                               3. W-W
% =======================================================================
FC_WW = WW(regionID_sortedByTracts, regionID_sortedByTracts);

vals_WW   = FC_WW(:);
num_pos   = sum(vals_WW > 0);
num_neg   = sum(vals_WW < 0);
num_zero  = sum(vals_WW == 0);
num_total = numel(vals_WW);

fprintf('W-W FC ：\n');
fprintf('  positive: %d (%.2f%%)\n',  num_pos,  100 * num_pos  / num_total);
fprintf('  negative: %d (%.2f%%)\n',  num_neg,  100 * num_neg  / num_total);
fprintf('  zero: %d (%.2f%%)\n\n', num_zero, 100 * num_zero / num_total);

figure('Position', [100, 100, 550, 500]);
imagesc(FC_WW);
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';
colormap(cmap);
caxis([-max(abs(FC_WW(:))), max(abs(FC_WW(:)))]);

hold on;
split_positions_WW = [10, 15, 34, 50];  
for i = 1:length(split_positions_WW)
    pos = split_positions_WW(i) + 0.5;

    line([0.5, 68.5], [pos, pos], 'Color', 'k', 'LineWidth', 0.5);

    line([pos, pos], [0.5, 68.5], 'Color', 'k', 'LineWidth', 0.5);
end
hold off;

xtick_positions_WW = [5, 13, 25, 43, 59];
ytick_positions_WW = [5, 13, 25, 43, 59];
xticks(xtick_positions_WW);
yticks(ytick_positions_WW);
xticklabels(tractsName);
yticklabels(tractsName);

ax = gca;
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
title('W-W FC');

save_name_WW = fullfile(savepath, 'PNC_WW_FCmean');
saveas(gcf, [save_name_WW '.fig']);
saveas(gcf, [save_name_WW '.jpg']);
print(gcf, [save_name_WW '.pdf'], '-dpdf', '-bestfit');