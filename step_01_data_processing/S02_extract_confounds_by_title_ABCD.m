function extract_confounds_by_title(filePath, savePath_bids, savePath_customCP, boldPath, dsegPath, brainmaskPath,titles)  
    % Read the TSV file  
    tbl = readtable(filePath, 'FileType', 'text', 'Delimiter', ' ');  

   %% get 6 motion parameters and their derivatives
    motion_params = [tbl.trans_x_mm, tbl.trans_y_mm, tbl.trans_z_mm, tbl.rot_x_degrees, tbl.rot_y_degrees, tbl.rot_z_degrees,...
         tbl.trans_x_mm_dt, tbl.trans_y_mm_dt, tbl.trans_z_mm_dt, tbl.rot_x_degrees_dt, tbl.rot_y_degrees_dt, tbl.rot_z_degrees_dt];
    motion_params_full = [tbl.trans_x_mm, tbl.trans_y_mm, tbl.trans_z_mm, tbl.rot_x_degrees, tbl.rot_y_degrees, tbl.rot_z_degrees,...
         tbl.trans_x_mm_dt, tbl.trans_y_mm_dt, tbl.trans_z_mm_dt, tbl.rot_x_degrees_dt, tbl.rot_y_degrees_dt, tbl.rot_z_degrees_dt, tbl.framewise_displacement];
    colNames_motion = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', ...
         'trans_x_derivative1','trans_y_derivative1','trans_z_derivative1', ...
         'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1', 'framewise_displacement'};
    % Convert array to table
    motionTable = array2table(motion_params_full);
    % Assign column names to the table
    motionTable.Properties.VariableNames = colNames_motion;
    % Compute power2 of motions parameters and dt 
    power2 = motion_params .^ 2;
    columnNames_motion_power2 = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2'...
         'trans_x_derivative1_power2','trans_y_derivative1_power2','trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'};
    % Convert array to table
    motionTable_power2 = array2table(power2);
    % Assign column names to the table
    motionTable_power2.Properties.VariableNames = columnNames_motion_power2;
        
    motionTable_merged = [motionTable, motionTable_power2];

   %% get csf mean signal
    bold_matrix = niftiread(boldPath);
    dseg_matrix = niftiread(dsegPath);
    % csf mask
    csf_mask = zeros(size(dseg_matrix));
    csf_mask(dseg_matrix==1) =1;
    % go over all timepoints (volumes)
    for i_ntpoints=1:size(bold_matrix,4)    
        curr_volume_data = bold_matrix(:,:,:,i_ntpoints);          
        curr_volume_data = curr_volume_data.*csf_mask; 
        bold_matrix_csf(:,:,:,i_ntpoints) = curr_volume_data;       
        clear curr_volume_data;
    end
    % get mean signal
    mean_csf_signal = zeros(size(bold_matrix_csf, 4), 1); 
    for i_t = 1:size(bold_matrix_csf, 4)
        mean_csf_signal(i_t) = mean(bold_matrix_csf(:, :, :, i_t), 'all');
    end

    %% get white matter mean signal
    % WM mask
    WM_mask = zeros(size(dseg_matrix));
    WM_mask(dseg_matrix==3) =1;
    % go over all timepoints (volumes)
    for i_ntpoints=1:size(bold_matrix,4)    % go over all timepoints (volumes)
        curr_volume_data = bold_matrix(:,:,:,i_ntpoints);          
        curr_volume_data = curr_volume_data.*WM_mask; 
        bold_matrix_WM(:,:,:,i_ntpoints) = curr_volume_data;       
        clear curr_volume_data;
    end    
    % get mean signal
    mean_WM_signal = zeros(size(bold_matrix_WM, 4), 1); 
    for i_t = 1:size(bold_matrix_WM, 4)
        mean_WM_signal(i_t) = mean(bold_matrix_WM(:, :, :, i_t), 'all');
    end

   %% get global mean signal
    % brain mask
    brain_matrix = niftiread(brainmaskPath);
    brain_mask = zeros(size(brain_matrix));
    brain_mask(brain_matrix==1) =1;
    % go over all timepoints (volumes)
    for i_ntpoints=1:size(bold_matrix,4)    
        curr_volume_data = bold_matrix(:,:,:,i_ntpoints);          
        curr_volume_data = curr_volume_data.*brain_mask; 
        bold_matrix_global(:,:,:,i_ntpoints) = curr_volume_data;       
        clear curr_volume_data;
    end
    % get mean signal
    mean_global_signal = zeros(size(bold_matrix_global, 4), 1); 
    for i_t = 1:size(bold_matrix_global, 4)
        mean_global_signal(i_t) = mean(bold_matrix_global(:, :, :, i_t), 'all');
    end

   %% get csf, WM, and global signal dt, power2, and dt's power2
    csf = mean_csf_signal; 
    csf2 = csf.^2; 
    csf_dt = diff(csf, 1, 1);
    csf_dt_nan = [zeros(1, size(csf, 2)); csf_dt];
    csf_dt_nan2 = csf_dt_nan.^2;

    WM = mean_WM_signal; 
    WM2 = WM.^2; 
    WM_dt = diff(WM, 1, 1);
    WM_dt_nan = [zeros(1, size(WM, 2)); WM_dt];
    WM_dt_nan2 = WM_dt_nan.^2;

    global_signal = mean_global_signal;
    global_signal2 = global_signal.^2;
    global_signal_dt = diff(global_signal, 1, 1);
    global_signal_dt_nan = [zeros(1, size(global_signal, 2)); global_signal_dt];
    global_signal_dt_nan2 = global_signal_dt_nan.^2;

   %% merge all confounds
    % Convert array to table
    colNames_brain = {'global_signal','white_matter','csf','global_signal_derivative1','white_matter_derivative1','csf_derivative1',...
        'global_signal_power2','white_matter_power2','csf_power2','global_signal_derivative1_power2','white_matter_derivative1_power2','csf_derivative1_power2'};
    dataTable_brain = array2table([global_signal, WM, csf, global_signal_dt_nan, WM_dt_nan, csf_dt_nan, global_signal2, WM2, csf2, global_signal_dt_nan2, WM_dt_nan2, csf_dt_nan2]);
    dataTable_brain.Properties.VariableNames = colNames_brain;

    dataTable_merged = [motionTable_merged, dataTable_brain];
    writetable(dataTable_merged, savePath_bids, 'FileType', 'text', 'Delimiter', '\t');  
    fprintf('new computed confounds files have been saved to: %s\n', savePath_bids); 
    
    %%
    % Initialize an empty table to store the extracted columns  
    extractedTbl = table;  
    % Iterate through all titles to extract the corresponding columns  
    for i = 1:length(titles)  
        % Find the column with a matching title  
        colIndex = find(strcmp(dataTable_merged.Properties.VariableNames, titles{i}));  
        % If a matching column is found, add it to the extracted table  
        if ~isempty(colIndex)  
            extractedTbl.(titles{i}) = dataTable_merged{:, colIndex};  
        else  
            warning('Column "%s" not found in the TSV file.', titles{i});  
        end  
    end  

    % Write the extracted table to a new TSV file  
    writetable(extractedTbl, savePath_customCP, 'FileType', 'text', 'Delimiter', '\t');  
    % Display the path of the output file  
    fprintf('Extracted columns have been saved to: %s\n', savePath_customCP);  

end
