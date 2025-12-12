%%
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
dataset_str = 'datasetsAveraged';

out_dir = [proj_dir '/Figure/Age_perdiction_feature_weight'];
% save_name = [out_dir '/' dataset_str '_sig_Haufe_FCmatrix_gw'];
save_name = [out_dir '/' dataset_str '_Haufe_FCmatrix_gw'];
%% get median feature weight for each FC
Results_Cell = g_ls([out_dir '/datasets/*_gw_Schaefer100.mat']);
% Results_Cell = g_ls([out_dir '/datasets_sig/*_Haufe_sigFCmatrix_gw.mat']);

for i_set = 1:length(Results_Cell)
  tmp = load(Results_Cell{i_set});
  FCmatrix_totalSets(:,:,i_set) = tmp.FCmatrix;
end

%%
FCmatrix_full = mean(FCmatrix_totalSets,3); 
save([save_name '.mat'], 'FCmatrix_full');

%% imagesec FCmatrix_full
load cmap0.mat
load Schaefer100_info.mat;
load JHU68_info_rename.mat
FCmatrix_rearranged = FCmatrix_full(regionID_sortedByNetwork, regionID_sortedByTracts);
figure('Position', [100, 100, 550, 500]);
imagesc(FCmatrix_rearranged); 
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';
colormap(cmap);
% caxis([-0.04, 0.04]);
caxis([-max(abs(FCmatrix_rearranged(:))), max(abs(FCmatrix_rearranged(:)))]);

    
%% add network lines
hold on;
col_separators = [10, 15, 34, 50]; 
for i = 1:length(col_separators)
    x = col_separators(i) + 0.5;
    line([x, x],[0.5, 100.5], 'Color', 'k', 'LineWidth', 0.5);  % Horizontal lines
end

row_separators = [17, 31, 46, 58, 63, 76]; 
for i = 1:length(row_separators)
    y = row_separators(i) + 0.5;
    line([0.5, 68.5], [y, y],  'Color', 'k', 'LineWidth', 0.5);  % Vertical lines
end
hold off;

%% add network legends
% white matter network
xtick_positions =  [5, 13, 25, 43, 59]; 
xticks(xtick_positions);
xtick_labels = tractsName;
xticklabels(xtick_labels);
% gray matter network
ytick_positions =  [8, 25, 39, 52,  61, 70, 90]; 
yticks(ytick_positions);
ytick_labels = networkName;
yticklabels(ytick_labels);
ax = gca;
ax.XTickLabelRotation = 90;
title('G-W feature weight matrix');
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; 

% save figure
saveas(gcf, [save_name '.fig']); 
saveas(gcf, [save_name '.jpg']); 
print(gcf, [save_name '.pdf'], '-dpdf', '-bestfit');

