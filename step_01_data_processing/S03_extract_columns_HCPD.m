function extract_columns_by_title_new2025(filePath, savePath_bids, savePath_customCP, titles)  
% Read the TSV file  
tbl = readtable(filePath, 'FileType', 'text', 'Delimiter', '\t');  

	if size(tbl,2)>10
		% Initialize an empty table to store the extracted columns  
		extractedTbl = table;  
		% Iterate through all titles to extract the corresponding columns  
		for i = 1:length(titles)  
			% Find the column with a matching title  
			colIndex = find(strcmp(tbl.Properties.VariableNames, titles{i}));  

			% If a matching column is found, add it to the extracted table  
			if ~isempty(colIndex)  
				extractedTbl.(titles{i}) = tbl{:, colIndex};  
			else  
				warning('Column "%s" not found in the TSV file.', titles{i});  
			end  
		end  
		% Write the extracted table to a new TSV file  
		writetable(extractedTbl, savePath_customCP, 'FileType', 'text', 'Delimiter', '\t');  
		% Display the path of the output file  
		fprintf('Extracted columns have been saved to: %s\n', savePath_customCP);  
	else
		%% Assume motion_params is a matrix where rows represent time points and columns represent different head motion parameters
		motion_params = [tbl.trans_x, tbl.trans_y, tbl.trans_z, tbl.rot_x, tbl.rot_y, tbl.rot_z, tbl.global_signal, tbl.white_matter, tbl.csf];
		
		%% 1. Compute the first-order difference (derivative1)
		derivative1 = diff(motion_params, 1, 1);
		% Create a matrix of zeros with the same size as the original data for padding
		derivative1_nan = [zeros(1, size(motion_params, 2)); derivative1];
		columnNames_derivative1_nan = {'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', 'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1'...
			 'global_signal_derivative1','white_matter_derivative1','csf_derivative1'};
		% Convert array to table
		dataTable_derivative1_nan = array2table(derivative1_nan);
		% Assign column names to the table
		dataTable_derivative1_nan.Properties.VariableNames = columnNames_derivative1_nan;

		%% 2. Compute power2 (squared head motion parameters)
		power2 = motion_params .^ 2;
		columnNames_power2 = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2'...
			 'global_signal_power2','white_matter_power2','csf_power2'};
		% Convert array to table
		dataTable_power2 = array2table(power2);
		% Assign column names to the table
		dataTable_power2.Properties.VariableNames = columnNames_power2;

		%% 3. Compute the first-order difference with power2 (derivative1_power2)
		diff_power2 = derivative1_nan.^ 2;
		columnNames_power2_diff_nan = {'trans_x_derivative1_power2', 'trans_y_derivative1_power2', 'trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'...
			 'global_signal_derivative1_power2','white_matter_derivative1_power2','csf_derivative1_power2'};
		% Convert array to table
		dataTable_diff_power2 = array2table(diff_power2);
		% Assign column names to the table
		dataTable_diff_power2.Properties.VariableNames = columnNames_power2_diff_nan;

		%% combine together
		fd = zeros(size(tbl,1),1);
		columnNames_fd = {'framewise_displacement'};
		dataTable_fd = array2table(fd);
		dataTable_fd.Properties.VariableNames = columnNames_fd;
			
		dataTable_merged = [tbl, dataTable_derivative1_nan, dataTable_power2, dataTable_diff_power2, dataTable_fd];
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
end