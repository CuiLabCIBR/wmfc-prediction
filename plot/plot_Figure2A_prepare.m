%%
clear;
ProjectFolder = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction';

targetStr_totalSets = {'HCPD', 'PNC', 'EFNY', 'CCNP'};
targetStr_totalSets = targetStr_totalSets';

targetDir_total = {'/HCPD/code/4th_prediction/s02_PLSprediction/nosmooth/age_withHaufe_Schaefer100',...
    '/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100',...
    '/EFNY/code/4th_prediction/s02_prediction/FC_nosmooth/age_withHaufe_Schaefer100_522',...
    '/CCNP/code/4th_prediction/nosmooth_n193_delete_under6/age_withHaufe_Schaefer100'};


dataCell_totalSets = cell(length(targetStr_totalSets)+1,6);
dataCell_totalSets(1,:) = {'tragetStr', 'R_gg_total', 'R_gw_total', 'partialR_gw_total', 'R_ww_total', 'partialR_ww_total'}; 
dataCell_totalSets(2:end,1) = targetStr_totalSets;

for i_set = 1:length(targetStr_totalSets)
    targetDir = targetDir_total{i_set};
    load([ProjectFolder, targetDir, '/partial_results_forBoxplot_5fold_sortedByCVtimes.mat']);
    dataCell_totalSets(i_set+1, 2:6) = [dataCell(2,2:end)];%% GG, GW, WW
end

save('age_dataCell_totalSets_forBoxplot.mat', 'dataCell_totalSets');