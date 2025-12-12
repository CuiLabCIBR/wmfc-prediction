%%
clear;
proj_dir =  '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';
addpath('/ibmgpfs/cuizaixu_lab/congjing/toolbox/PANDA_1.3.0_64');
addpath([proj_dir '/post_prediction/Ageprediction_feature_weight']);

datasets_str = {'HCPD', 'PNC', 'EFNY', 'CCNP'};

datasets_dir = { ...
    [proj_dir '/HCPD/code/4th_prediction/s02_PLSprediction/nosmooth/age_withHaufe_Schaefer100/HCPD_interview_age'], ... % HCPD
    [proj_dir '/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100/PNC_age'], ... % PNC
    [proj_dir '/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth/age_withHaufe_Schaefer100_522/brainproject522_age'], ... % EFNY
    [proj_dir '/CCNP/code/4th_prediction/nosmooth_n193_delete_under6/age_withHaufe_Schaefer100/CCNP_age'] ... % CCNP
};


out_dir = [proj_dir '/Figure/Age_perdiction_feature_weight/datasets'];
%%
for i_set = 1:length(datasets_str)
    current_set = datasets_str{i_set}
    data_dir = datasets_dir{i_set};
    save_name = [out_dir '/' current_set '_Haufe_FCmatrix_gg_Schaefer100'];
    
   %% get median weight for each FC with 101CV
    Results_Cell = g_ls([data_dir '/RegressCovariates_RandomCV/Time_*/GGFC/Fold_*_Score.mat']);
    clear w_Brain_Overall_101CV
    for i_cv = 1:length(Results_Cell)
      tmp = load(Results_Cell{i_cv});
      w_Brain_Overall_101CV(i_cv, :) = tmp.w_Brain_Haufe;
    end
    w_Brain_median_101CV = median(w_Brain_Overall_101CV, 1);
    
  
   %% 
    % w_Brain_median_101CV = median(w_Brain_Overall_101CV);
    VisualizeFolder = [data_dir '/WeightVisualize_Overall_RandomCV'];
    mkdir(VisualizeFolder);
    save([VisualizeFolder '/w_Brain_median_101CV_GG_Haufe.mat'], 'w_Brain_median_101CV');

   %% convert from FC vector to FC matrix
    lower_tri_vec = w_Brain_median_101CV; 
    m=100; n=100;
    % Initialize an empty m x n matrix
    FCmatrix = zeros(m, n);

    if strcmpi(current_set, 'EFNY')
        % ===== EFNY：Python ( (2,1),(3,1),(3,2),... ) =====
        [ri, ci] = find(tril(true(m, n), -1));       
        [~, ord] = sortrows([ri, ci], [1 2]);        
        lin_idx = sub2ind([m n], ri(ord), ci(ord));  
        FCmatrix(lin_idx) = lower_tri_vec(:);
    else
        % ===== HCPD/PNC/CCNP：MATLAB  =====
        tril_idx = find(tril(true(m, n), -1));       
        FCmatrix(tril_idx) = lower_tri_vec(:);
    end
    
    FCmatrix_full = FCmatrix + FCmatrix';
    % FCmatrix_full(logical(eye(m))) = 0;
    save([save_name '.mat'], 'FCmatrix', 'FCmatrix_full');
    
    %% sorted FCmatrix by networks and imagesec 
    load Schaefer100_info.mat;
    FCmatrix_rearranged = FCmatrix_full(regionID_sortedByNetwork, regionID_sortedByNetwork);
    % set upper triangle with 0
    % FCmatrix_rearranged(triu(true(size(FCmatrix_rearranged)))) = 0;
    figure('Position', [100, 100, 550, 500]);
    imagesc(FCmatrix_rearranged); 
    colorbar
    caxis([-max(abs(FCmatrix_rearranged(:))), max(abs(FCmatrix_rearranged(:)))]);
    colormap(redblue);
    
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
    title(['G-G weight matrix of ' current_set]);
    ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; 
    % save figure
    saveas(gcf, [save_name '.fig']); 
    saveas(gcf, [save_name '.jpg']); 

end
