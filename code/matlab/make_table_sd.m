% Define the root directory (can be the parent directory containing all folders)
function make_table_sd(specif)

rootDir = '/Volumes/Newsmy 1/phd/cepstral_distance/results/matlab/simulation/mean_sd';

% Get all folders containing "Martin" in the root directory

spec = strcat('*',specif,'*');
folders = dir(fullfile(rootDir, spec));

% Filter to keep only folders
folders = folders([folders.isdir]);

folderNames = {folders.name};

% The number of folders
K = size(folderNames, 2);

% Define folder path
table = [];
for i = 1:K
folderPath = folderNames{i};

% Get the number after the first 'N' in the folder name
folderName = folderPath; % Folder name
N_index = strfind(folderName, 'N'); % Find the position of 'N'
dash_index = strfind(folderName, 'SNR'); % Find the position of 'SNR'
number_after_N = str2double(folderName(N_index+1:dash_index-1)); % Extract numbers

% Extract the part after '-1'
dash_index_2 = strfind(folderName, specif);
a = length(specif);
% Find the position of '-1'
suffix_after_dash = folderName(dash_index_2+a:end); % Extract the part after '-1'

% Initialize an empty array to hold all ResCepNulling data
allResCepNulling = [];

% Define the file names to be merged
fileNames = {'0.01resultsFDR.mat', '0.05resultsFDR.mat', 'resultsCN.mat', 'resultsWP.mat'};

% Iterate through each file and load the ResCepNulling data
for i = 1:length(fileNames)
    data = load(fullfile(folderPath, fileNames{i}), 'ResCepNulling'); % 加载指定变量
    
    allResCepNulling = [allResCepNulling; data.ResCepNulling];
end

% Add a column in allResCepNulling, the value of this column is number_after_N
numberColumn = num2cell(number_after_N * ones(size(allResCepNulling, 1), 1));
% Create a cell array column with the same number of rows as myCellArray, where each item is 'Martin'
newColumn = cell(size(allResCepNulling, 1), 1);
newColumn(:) = {suffix_after_dash};
allResCepNulling = [allResCepNulling, numberColumn, newColumn];
table = [table; allResCepNulling];

end


% Convert the cell array to a table for sorting
myTable = cell2table(table, 'VariableNames', {'Method', 'Mean', 'Variance', 'N', 'Weight'});

% Sort by the fourth and fifth columns in descending order
sortedTable = sortrows(myTable, {'Weight', 'N'}, {'descend', 'descend'});


% Save the table as an Excel file
name = strcat('simulation_mean_sd_',specif,'.xlsx');
writetable(sortedTable, name);

end
