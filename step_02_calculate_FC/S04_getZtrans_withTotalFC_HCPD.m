%% 
clear;
load subid_final_hcpd.mat;
subList = subid_final_HCPD;
datapath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/results';
outpath = '/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/HCPD/results';
savepath = [outpath '/FC_total_nosmooth'];

GM_network_num = 100; %% Yeo7
WM_network_num = 68; %% ICBM

total_FCmatrix_WW = zeros(WM_network_num,WM_network_num,length(subList));
total_FCmatrix_GW = zeros(GM_network_num,WM_network_num,length(subList));
total_FCmatrix_GG = zeros(GM_network_num,GM_network_num,length(subList));

%%
for i_sub=1:length(subList)
    i_sub
    sub_current = subList{i_sub};
    %% 1. WM-WM
    load([datapath '/FCmatrix_individual_Schaefer100_nosmooth/sub-' sub_current '/sub-' sub_current '_rICBM_60p_WW_FC.mat']);
   % load([datapath '/FCmatrix_individual_Schaefer100/' sub_current '/' sub_current '_rICBM_60p_WW_FC.mat']);
    WW_FC = WW_FCmatrix;
    %%% fisher z
    WW_FC_z = atanh(WW_FC);
    FCmatrix = WW_FC_z;
    [m, n] = size(FCmatrix);  
    %%% convert matrix to vector
    FCvector = tril(FCmatrix,-1);
    FCvector_WW = FCvector(tril(true(m), -1));  
    total_FCvector_WW(i_sub, :) = FCvector_WW;
    total_FCmatrix_WW(:,:,i_sub) = FCmatrix;
       
    %% 2. GM-WM
    load([datapath '/FCmatrix_individual_Schaefer100_nosmooth/sub-' sub_current '/sub-' sub_current '_Schaefer100_rICBM_60p_GW_FC.mat']);
%    load([datapath '/FCmatrix_individual_Schaefer100/' sub_current '/' sub_current '_Schaefer100_rICBM_60p_GW_FC.mat']);
    GW_FC = GW_FCmatrix;
    %%% fisher z
    GW_FC_z = atanh(GW_FC);
    FCmatrix = GW_FC_z;
    %%% convert matrix to vector
    FCvector = reshape(FCmatrix, 1, []);
    FCvector_GW = FCvector;
    total_FCvector_GW(i_sub, :) = FCvector_GW;
    total_FCmatrix_GW(:,:,i_sub) = FCmatrix;
    
   %% 3. GM-GM
    load([datapath '/FCmatrix_individual_Schaefer100_nosmooth/sub-' sub_current '/sub-' sub_current '_Schaefer100_GG_FC.mat']);
%    load([datapath '/FCmatrix_individual_Schaefer100/' sub_current '/' sub_current '_Schaefer100_GG_FC.mat']);
    GG_FC = GG_FCmatrix;
    %%% fisher z
    GG_FC_z = atanh(GG_FC);
    FCmatrix = GG_FC_z;
    [m, n] = size(FCmatrix);  
    %%% convert matrix to vector
    FCvector = tril(FCmatrix,-1);
    FCvector_GG = FCvector(tril(true(m), -1)); 
    total_FCvector_GG(i_sub, :) = FCvector_GG;
    total_FCmatrix_GG(:,:,i_sub) = FCmatrix;
    
end

save([savepath '/total_FCmatrix_GW.mat'], 'total_FCmatrix_GW');
save([savepath '/total_FCvector_GW.mat'], 'total_FCvector_GW');
save([savepath '/total_FCmatrix_GG.mat'], 'total_FCmatrix_GG');
save([savepath '/total_FCvector_GG.mat'], 'total_FCvector_GG');
save([savepath '/total_FCmatrix_WW.mat'], 'total_FCmatrix_WW');
save([savepath '/total_FCvector_WW.mat'], 'total_FCvector_WW');

%%% remove NAN value in the matrix
A = total_FCvector_GW;
A(isnan(A)) = 0;
writematrix(A, [savepath '/total_FCvector_GW.txt']);
A = [];
A = total_FCvector_GG;
A(isnan(A)) = 0;
writematrix(A, [savepath '/total_FCvector_GG.txt']);
A = [];
A = total_FCvector_WW;
A(isnan(A)) = 0;
writematrix(A, [savepath '/total_FCvector_WW.txt']);

writecell(cellstr(subList), fullfile(savepath,'subjects_order.txt'));