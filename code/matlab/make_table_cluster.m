function make_table_cluster(specif)

% Get the paths of folders containing "SNR-1" in their names
folders = '/Volumes/Newsmy 1/phd/cepstral_distance/code/matlab/Simulation';

% Get all folders containing "Martin" in the root directory

spec = strcat('*',specif,'*');
folders = dir(fullfile(folders, spec));

% Define the order of the Method column and repeat it twice
method_order = {'rect', 'hann', 'log-Periodogram', 'BIC', 'KSF', 'MRI', 'FDR_0.01', 'FER_0.01', 'FDR_0.05', 'FER_0.05'};
method_repeated = repmat(method_order, 1, 2);

% Initialize the final table to store results
final_data = [];

for i = 1:length(folders)
    folder_name = folders(i).name;
    mat_file_path = fullfile(folder_name, 'results.mat');
    
    if exist(mat_file_path, 'file')
        % Load data from the results.mat file
        data = load(mat_file_path, 'sim_ID', 'sim_martin');
        
        % Check if both sim_ID and sim_martin variables exist
        if isfield(data, 'sim_ID') && isfield(data, 'sim_martin')
            % Extract sim_ID and sim_martin data
            sim_ID = data.sim_ID;
            sim_martin = data.sim_martin;
            
            % Combine data and add the Weight column
            sim_ID_data = [num2cell(sim_ID), repmat({'Identity'}, size(sim_ID, 1), 1)];
            sim_martin_data = [num2cell(sim_martin), repmat({'Martin'}, size(sim_martin, 1), 1)];
            
            % Merge sim_ID and sim_martin data
            combined_data = [sim_ID_data; sim_martin_data];
            
            % Add the Method column
            num_rows = size(combined_data, 1);
            combined_data = [combined_data, method_repeated(1:num_rows)'];
            
            % Extract the numeric value for N from the folder name
            n_value = regexp(folder_name, '(?<=N)(\d+)(?=SNR)', 'match', 'once');
            n_column = repmat({str2double(n_value)}, num_rows, 1);
            
            % Combine all data into a table
            combined_table = table(combined_data(:, 1), combined_data(:, 2), combined_data(:, 3), n_column, ...
                'VariableNames', {'sim_metrics', 'Weight', 'Method', 'N'});
            
            % Append the table to the final dataset
            final_data = [final_data; combined_table];
        else
            % Show a warning if the required variables are missing
            warning('Missing variables in %s. Skipping...', mat_file_path);
        end
    else
        % Show a warning if the file is not found
        warning('File not found: %s', mat_file_path);
    end
end

% Save the final data to a MAT file or export it as a CSV

name = strcat(specif,'_similarity_index.xlsx');
writetable(final_data, name);
end


