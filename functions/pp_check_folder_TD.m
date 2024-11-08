function pp_check_folder_TD(varargin)
% pp_check_folder_TD checks if the specified folder paths exist.
% If a path does not exist, the function prompts the user to create it.
%
% Usage:
%   pp_check_folder_TD('Path1', 'Path2', ...)
%
% Inputs:
%   varargin - any number of folder paths as strings
% Tim Dressler, 07.11.2024

for i = 1:numel(varargin)
    folderPath = varargin{i};
    if ~exist(folderPath)
        fprintf('The folder "%s" does not exist.\n', folderPath);
        %ask whether to create the folder
        answer = input('Would you like to create this folder? (y/n): ', 's');
        if lower(answer) == 'y'
            try
                mkdir(folderPath);
                fprintf('Folder "%s" was successfully created\n', folderPath);
            catch
                fprintf('Error: Could not create folder "%s"\n', folderPath);
            end
        else
            fprintf('Folder "%s" was not created\n', folderPath);
        end
    else
        fprintf('Folder "%s" exists\n', folderPath);
    end
end
end
