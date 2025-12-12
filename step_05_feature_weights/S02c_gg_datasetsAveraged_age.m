%%
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
dataset_str = 'datasetsAveraged';

out_dir = [proj_dir '/Figure/Age_perdiction_feature_weight'];
save_name = [out_dir '/' dataset_str '_Haufe_FCmatrix_gg'];
% save_name = [out_dir '/' dataset_str '_sig_Haufe_FCmatrix_gg'];

%% get median feature weight for each FC
Results_Cell = g_ls([out_dir '/datasets/*_gg_Schaefer100.mat']);
% Results_Cell = g_ls([out_dir '/datasets_sig/*_Haufe_sigFCmatrix_gg.mat']);

for i_set = 1:length(Results_Cell)
  tmp = load(Results_Cell{i_set});
  FCmatrix_totalSets(:,:,i_set) = tmp.FCmatrix_full;
end

%%
FCmatrix_full = mean(FCmatrix_totalSets,3);
save([save_name '.mat'], 'FCmatrix_full');

%% imagesec FCmatrix_full
load cmap0.mat
load Schaefer100_info.mat;
FCmatrix_rearranged = FCmatrix_full(regionID_sortedByNetwork, regionID_sortedByNetwork);

figure('Position', [100, 100, 550, 500]);
imagesc(FCmatrix_rearranged); 
cb = colorbar;
cb.FontSize = 14;
cb.FontName = 'Arial';
colormap(cmap);
caxis([-max(abs(FCmatrix_rearranged(:))), max(abs(FCmatrix_rearranged(:)))]);

    
%% add network lines
hold on;
% Define the positions (between rows/columns) where separator lines are needed
split_positions = [17, 31, 46, 58, 63, 76]; 
% Draw horizontal lines between specified rows
for i = 1:length(split_positions)
    y = split_positions(i) + 0.5;
    line([0.5, 100.5], [y, y], 'Color', 'k', 'LineWidth', 0.5);
end
% Draw vertical lines between specified columns
for i = 1:length(split_positions)
    x = split_positions(i) + 0.5;
    line([x, x], [0.5, 100.5], 'Color', 'k', 'LineWidth', 0.5);
end
hold off;

%% add network legends
xtick_positions =  [8, 25, 39, 52,  61, 70, 90]; 
xticks(xtick_positions);
xtick_labels = networkName;
xticklabels(xtick_labels);
ytick_positions =  [8, 25, 39, 52, 61, 70, 90]; 
yticks(ytick_positions);
ytick_labels = networkName;
yticklabels(ytick_labels);
ax = gca;
ax.XTickLabelRotation = 90;
title('G-G feature weight matrix');
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; 

% save figure
saveas(gcf, [save_name '.fig']); 
saveas(gcf, [save_name '.jpg']); 
print(gcf, [save_name '.pdf'], '-dpdf', '-bestfit');

