function pp_clean_up_folder_TD(folderPath)
    % pp_clean_up_folder_TD - Moves files to a timestamped archive subfolder within an "archive" folder
    %
    %
    % Usage:
    %   pp_clean_up_folder_TD('Path')
    %
    % Input:
    %   folderPath - The path of the folder to check and clean up.
    %
    % Description:
    %   This function checks if a specified folder is empty. If the folder
    %   contains any files, it ensures an "archive" folder exists within the
    %   specified folder, creates a new subfolder named archive_DD_MM_YYYY_TIME
    %   inside the archive folder, and moves all files from the original folder 
    %   to this new archive subfolder.
    %
    % Tim Dressler, 12.11.24

    %check if the folder exists
    if ~isfolder(folderPath)
        error('Specified folder does not exist');
    end

    %get list of files in the folder (excluding '.' and '..')
    files = dir(folderPath);
    files = files(~ismember({files.name}, {'.', '..', 'archive'}));

    %if folder is not empty, proceed with archiving
    if ~isempty(files)
        %define the archive folder path
        archiveFolderPath = fullfile(folderPath, 'archive');
        
        %create the archive folder if it does not exist
        if ~exist(archiveFolderPath, 'dir')
            mkdir(archiveFolderPath);
        end

        %get current date and time in DD_MM_YYYY_TIME format
        currentDateTime = datestr(now, 'dd_mm_yyyy_HH-MM-SS');

        %create the timestamped archive subfolder within the archive folder
        timestampedArchiveFolder = fullfile(archiveFolderPath, ['archive_' currentDateTime]);
        
        %create the timestamped archive folder
        mkdir(timestampedArchiveFolder);

        %move files to the timestamped archive folder
        for k = 1:length(files)
            sourceFile = fullfile(folderPath, files(k).name);
            destinationFile = fullfile(timestampedArchiveFolder, files(k).name);

            %move the file
            movefile(sourceFile, destinationFile);
        end

        fprintf('Files moved to: %s\n', timestampedArchiveFolder);
    else
        fprintf('The folder is empty. No action taken.\n');
    end
end
