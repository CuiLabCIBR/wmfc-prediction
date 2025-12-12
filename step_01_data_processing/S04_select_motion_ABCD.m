%% at least 2 run with: 
% 1) mean FD<=0.5mm; 
% 2) more than 40%TR with FD<=0.2mm; 
% 3) timepoints>=half original time-series length
clear;
bids_path = '/ibmgpfs/cuizaixu_lab/zhaoshaoling/HCP_WMfMRI/ABCD/xcpd0.7.1rc5/bids';

load subID_QSIprep.mat; 
load subID_QCpassed_inABCC.mat;
subList = intersect(subID_QSIprep, subID_QCpassed_inABCC);

nSub = length(subList)
nMaxRuns = 2;
valid_subID = [];

% T_name = 'ABCD_MotionInfo_bids.csv';
% validT_name = 'ABCD_MotionInfo_bids_QCpassed.csv';
T_name = 'ABCD_MotionInfo_n5959_corrected_withTPchecked.csv'; %%TP: timepoints
validT_name = 'ABCD_MotionInfo_n5959_FDpassed_by_fmriprep_corrected_withTPchecked.csv';

%% each run has 4 cols: meanFD, FD02, timepoint, valid index
MotionInfo_totalSubj = nan(length(subList), 4*nMaxRuns); %% max run number is 4
subList_QCpassed = cell(length(subList), 1);

for i_sub=1:length(subList)
% for i_sub=1:50
    i_sub
    sub_current = subList{i_sub};
    dir_current = [bids_path '/sub-' sub_current '/func'];
   %%% check if motion file is existed
    MotionInfo_currentSubject = [];
%     motion_flag = dir([dir_current '/' sub_current '_ses-01_task-rest_*_desc-confounds_timeseries.tsv']);
    motion_flag = dir([dir_current '/sub-' sub_current '*_desc-includingFD_motion.tsv']);
	motion_exist = ~isempty(motion_flag);
    if motion_exist==1
        for i_run = 1:length({motion_flag.name})
            if i_run > 2
                break;
            end
            motion_current_run = readtable([motion_flag(i_run).folder '/' motion_flag(i_run).name], 'FileType', 'text', 'Delimiter', ' ');
%             FD = cellfun(@str2double, motion_current_run.framewise_displacement);
            FD = motion_current_run .framewise_displacement;
%             FD(isnan(FD)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            FD_valid = FD(2:end); %% exclude the first timepoint
            meanFD = nanmean(FD_valid);
            FD02=length(find(FD_valid<=0.2))/length(FD_valid); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%             if meanFD<=0.5&&FD02>=0.4
            if meanFD<=0.5&&FD02>=0.4&&length(FD)>=383/2
                subList_QCvalid_Index = 1;
            else 
                subList_QCvalid_Index = 0;
            end
            MotionInfo_currentSubject(1, (i_run*4-3)) = meanFD;
            MotionInfo_currentSubject(1, (i_run*4-2)) = FD02;
            MotionInfo_currentSubject(1, (i_run*4-1)) = length(FD);
            MotionInfo_currentSubject(1, (i_run*4)) = subList_QCvalid_Index;
        end
        FD_total_run = MotionInfo_currentSubject(1,1:4:end); 
        meanFD_currentSubj = mean(FD_total_run);
        FD02_total_run = MotionInfo_currentSubject(1,2:4:end); 
        FD02_currentSubj = mean(FD02_total_run);
    else 
        meanFD_currentSubj = 999;
        FD02_currentSubj = 0;
    end
    meanMotion_totalSubj(i_sub,1) = meanFD_currentSubj;
    meanMotion_totalSubj(i_sub,2) = FD02_currentSubj;
    MotionInfo_totalSubj(i_sub,1:size(MotionInfo_currentSubject,2)) = MotionInfo_currentSubject;
end

%% write motion info into a table, each run has 4 cols: meanFD, FD02, timepoint, valid index
T = [cell2table(subList), array2table(MotionInfo_totalSubj), array2table(meanMotion_totalSubj)];
% Number of run
A = MotionInfo_totalSubj;
% Number of runs (each run = 4 columns)
nRun = size(A,2) / 4;  
% Generate variable names
varNames = {'subList'};  % Start with sublist column
for r = 1:nRun
    varNames{end+1} = sprintf('meanFD_run%d', r);
    varNames{end+1} = sprintf('FD02ratio_run%d', r);
    varNames{end+1} = sprintf('timepoints_run%d', r);
    varNames{end+1} = sprintf('QCvalid_run%d', r);
end
varNamesT = [varNames, {'meanFD', 'FDratio_lessthan_02'}];
T.Properties.VariableNames = varNamesT ;
writetable(T, T_name);

%% QC pass for total run
QCvalid_total_run = nansum(MotionInfo_totalSubj(:, [4,8]), 2); %% each run has 4 cols: meanFD, FD02, timepoint, valid index
valid_id = find(QCvalid_total_run>=2);
valid_subID = subList(valid_id);
validT = T(valid_id, :);
writetable(validT, validT_name);

