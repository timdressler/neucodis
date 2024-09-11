% pp_2_pre_proc.m
%
% Preproccessing pipeline for eeg-data 2/3. 
%CHECK: Complete description after adaption
%
% Originally by Dr. Cassie Short & Prof. Dr. Stefan Debener
% Adapted from Tim Dressler, 11.09.2024

mainpath = 'C:\Users\short\Documents\ERPSamplingStudy'; % main path
path_eeglab = [mainpath, '\EEGLAB']; % where eeglab is located
path_chanloc = [mainpath, '\RecordingInfo']; % where chanloc struct is located
path_data = [mainpath, '\PreprocessedData\Chanlocsadded']; % where channellocationadded data is located
path_preprocessed = [mainpath, '\PreprocessedData']; % where pre-processed data is saved
path_data = [path_preprocessed, '\primary'];
path_condspecific = [mainpath, '\ConditionSpecificData']; % where condition specific datasets are saved
path_StdDevImages = [mainpath, '\StdDevImages']; % where standard deviation maps are saved for each participant, used to identify bad channels
addpath('C:\Users\short\Documents\ERPSamplingStudy\Scripts');

bsl_fork = {-200; -100};
ref_fork = {'Mas'; 'Avg'; 'CSD'};

responses = {'S  1', 'S  2', 'S  3', 'S  4', 'S  5', 'S  6', 'S  7'};
conds = {'fear', 'anger', 'disgust', 'happiness', 'neutral', 'sadness', 'surprise'};

win_fork = {500,200; 500,300; 600,200; 600,300; 600,600; 700,200; 700,300;...
    700,600; 450,100; 'GAV', 200; 'SAV', 200};
elec_fork = {{'CP1', 'CP2', 'Pz', 'P3', 'P4'}; {'P3', 'P4', 'CP1', 'CP2'}; {'P3','Pz','P4'};...
    {'Fz','Cz', 'Pz'}; {'CP1', 'CP2'}; {'Cz'}; {'Pz'};{'around_peak'}};

midline = {'Fpz', 'Fz', 'Cz', 'Pz', 'Oz', 'Iz'};

storepath = 'C:\Users\short\Documents\ERPSamplingStudy\Analysis\'; % main path

% Get a list of all subfolders in the 'primary' directory
subfolderStruct = dir(fullfile(path_data, '*'));  % Use full path and wildcard for subfolders
subfolders = {subfolderStruct([subfolderStruct.isdir]).name};  % Filter out non-folders and extract names

files = {};
subs = {};
for k = 1:length(subfolders)
    folderName = subfolders{k};  % Correctly extract the folder name as a string
    % Skip '.' and '..' directories to be ignored in the folder
    if strcmp(folderName, '.') || strcmp(folderName, '..')
        continue;
    end
    % Full path to the subfolder
    subfolderPath = fullfile(path_preprocessed, 'primary', folderName);
    % List of .set files in the subfolder
    setFiles = dir(fullfile(subfolderPath, '*_postICA_interpolated.set'));
    % Append these files (with their paths) to the 'files' array
    for j = 1:length(setFiles)
        files{end + 1} = fullfile(subfolderPath, setFiles(j).name);
        subs{end + 1} = folderName;
    end
end

SAVpeaks = {}; % necessary to store them seperately per subject at this step

for i = 1:length(files)

    cd(path_eeglab);
    close all; clc;
    eeglab; % start eeglab (and restart it after every iteration of the loop)

    EEG = pop_loadset('filename',files{i});

    % epoch data around events
    events = {'S 18', 'S 19', 'S 28', 'S 29', 'S 38', 'S 39', 'S 48', 'S 49', 'S 56', 'S 58', 'S 68', 'S 69', 'S 78', 'S 79'};
    EEG = pop_epoch(EEG, events, [-0.2 0.8]);

    rejsublist = {}; % to be filled with participant numbers of those who were not considered for the analysis
    rejsub = 0; % if a pp is rejected, this will be set to 1 later and the data will not be saved

    EEG_ori = EEG;

    for b = 1:length(bsl_fork)
        timerange = [bsl_fork{b} 0]; % baseline time range
        EEG = pop_rmbase(EEG_ori, timerange); % baseline correction

        %remove bad epochs
        datachan = [1:size(EEG.data,1)];
        list_props = epoch_properties(EEG,datachan);
        marked_trials = find(min_z(list_props));
        bad_trials = zeros(1,EEG.trials);
        bad_trials(marked_trials) = 1; % index to bad trials in the data
        EEG = pop_rejepoch(EEG,bad_trials,0); % bad trials are removed
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');

        for cn = 1:length(conds)
            ntrials = sum(strcmp({EEG.event.type},events{2*cn-1})) + ...
                sum(strcmp({EEG.event.type},events{2*cn}));
            if ntrials < 114*0.25
                rejsub = 1; % if a participant has less than 25% of viable trials for a condition, their data will not be included further
            end
        end

        if rejsub == 0 % = online execute this if participant is not rejected
            for r = 1:length(ref_fork) % reference fork
                if strcmp(ref_fork{r},'Mas') % mastoids reference
                    EEG = pop_reref(ALLEEG(2), [39 40],'refloc',struct('labels',{'Cz'},'ref',{'Cz'},'theta',{0},'radius',{0},'X',{6.1232e-17},'Y',{0},'Z',{1},'sph_theta',{0},'sph_phi',{90},'sph_radius',{1},'type',{''},'urchan',{[]},'datachan',{0}));
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'overwrite','on','gui','off');
                    crm = {EEG.chanlocs.labels}; % channels reference mastoids; this will be relevant later
                elseif strcmp(ref_fork{r}, 'Avg') % average reference
                    EEG = pop_reref(ALLEEG(2), [],'refloc',struct('labels',{'Cz'},'ref',{'Cz'},'theta',{0},'radius',{0},'X',{6.1232e-17},'Y',{0},'Z',{1},'sph_theta',{0},'sph_phi',{90},'sph_radius',{1},'type',{''},'urchan',{[]},'datachan',{0}));
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'overwrite','on','gui','off');
                    cra = {EEG.chanlocs.labels}; % channels reference average; this will be relevant later
                elseif strcmp(ref_fork{r},'CSD') % common source density reference
                    EEG = pop_reref(ALLEEG(2), [],'refloc',struct('labels',{'Cz'},'ref',{'Cz'},'theta',{0},'radius',{0},'X',{6.1232e-17},'Y',{0},'Z',{1},'sph_theta',{0},'sph_phi',{90},'sph_radius',{1},'type',{''},'urchan',{[]},'datachan',{0})); %reref so Cz is back in the data fo the CSD
                    % Find the index of A1 and A2
                    a1Index = find(strcmp({EEG.chanlocs.labels}, 'A1'));
                    a2Index = find(strcmp({EEG.chanlocs.labels}, 'A2'));
                    % Mirror A2 coordinates for A1
                    EEG.chanlocs(a1Index).X = -EEG.chanlocs(a2Index).X; % Invert X-coordinate
                    EEG.chanlocs(a1Index).Y = EEG.chanlocs(a2Index).Y;  % Keep Y-coordinate the same
                    EEG.chanlocs(a1Index).Z = EEG.chanlocs(a2Index).Z;  % Keep Z-coordinate the same
                    [nChannels, nTimePoints, nEpochs] = size(EEG.data);
                    csdData = zeros(nChannels, nTimePoints, nEpochs);  % Initialize container for CSD data
                    for epoch = 1:nEpochs
                        % Extract data for the current epoch
                        epochData = EEG.data(:, :, epoch);
                        % Reshape epoch data to 2D
                        epochData2D = reshape(epochData, nChannels, nTimePoints);
                        % Apply CSD transformation
                        csdEpochData = laplacian_perrinX(epochData2D, [EEG.chanlocs.X], [EEG.chanlocs.Y], [EEG.chanlocs.Z]);
                        % Store the transformed data
                        csdData(:, :, epoch) = csdEpochData;
                    end
                    % Update EEG.data with the CSD-transformed data
                    EEG.data = csdData;
                    % Update the EEG structure in ALLEEG
                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'overwrite', 'on', 'gui', 'off');
                    crc = {EEG.chanlocs.labels}; % Channels reference CSD; this will be relevant later
                end

                ERPc = [];
                for c = 1:length(conds) % conditions
                    EEG = pop_selectevent(ALLEEG(3), 'type',{events{2*c-1} events{2*c}},'deleteevents','off','deleteepochs','on','invertepochs','off');
                    EEG = pop_eegfiltnew(EEG, [],30,220,0,[],0);
                    if ~exist([storepath 'fork1-2'], 'dir')
                        mkdir([storepath 'fork1-2'])
                    end
                    currentSub = subs{i}; % Store the current participant in a variable
                    bsl_num = bsl_fork{b}; % Get the current bsl_fork value

                    if isequal(bsl_num, -200) && strcmp(ref_fork{r}, 'Mas')
                        % Save the outcome for the current participant (bsl_fork option 1 with ref_fork option 1)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    elseif isequal(bsl_num, -200) && strcmp(ref_fork{r}, 'Avg')
                        % Save the outcome for the current participant (bsl_fork option 1 with ref_fork option 2)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    elseif isequal(bsl_num, -200) && strcmp(ref_fork{r}, 'CSD')
                        % Save the outcome for the current participant (bsl_fork option 1 with ref_fork option 3)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    elseif isequal(bsl_num, -100) && strcmp(ref_fork{r}, 'Mas')
                        % Save the outcome for the current participant (bsl_fork option 2 with ref_fork option 1)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    elseif isequal(bsl_num, -100) && strcmp(ref_fork{r}, 'Avg')
                        % Save the outcome for the current participant (bsl_fork option 2 with ref_fork option 2)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    elseif isequal(bsl_num, -100) && strcmp(ref_fork{r}, 'CSD')
                        % Save the outcome for the current participant (bsl_fork option 2 with ref_fork option 3)
                        outcomeFile = sprintf('%s%s_%s_b%d%sr%s.set', storepath, 'fork1-2', currentSub, bsl_num, conds{c}, ref_fork{r});
                    end

                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'overwrite', 'on', 'savenew', outcomeFile, 'gui', 'off');
                    ERPc(c,:,:) = mean(EEG.data,3);
                end

                % The timing fork will take place only at the LPP parametrisation, but we
                % already need notes of peak locations for subject average and grand
                % average for later
                GAVlist{b,r}(i,:,:,:) = ERPc; % creates b x r (nBaselines, nReferences) cells that contain ERPs per condition and participant; used in electrode: GAV peak path later
                SAV = squeeze(mean(ERPc,1)); % saves participant i's average ERP across conditions
                Dpoints_sav = round(eeg_lat2point([0.3 1], [1 1], EEG.srate, [EEG.xmin EEG.xmax]));

                max_amp = max(max(SAV(:,Dpoints_sav(1):Dpoints_sav(2)))); % find the max amplitude within the time window of interest
                [chan lat_max_amp] = find(SAV == max_amp); % find the indicies for that amplitude
                lat_max_amp_ms = EEG.times(lat_max_amp); % convert the latency from sample point to ms
                SAVpeaks(end+1,:) = {[currentSub,'_b',num2str(bsl_fork{b}),'r',ref_fork{r}],lat_max_amp_ms,chan}; % saves {participantID & path, timing of peak, channelindex of peak)
            end
        else
            rejsublist{end+1} = currentSub;
        end
    end

    idx = strfind(files{i}, '_postICA_interpolated.set') - 1;
    newfilename = [files{i}(1:idx), '_postICA_epoched_pruned.set'];

    EEG = pop_saveset(EEG, 'filename', newfilename);
end

%% GAV peaks
%
GAVpeaks = {'b-100rMas';'b-200rMas';'b-100rAvg';'b-200rAvg';'b-100rCSD';'b-200rCSD'}; % hardcoded for simplicity since b x r is only 2x2
for g = 1:numel(GAVlist)
    GAV = squeeze(squeeze(mean(mean(GAVlist{g},2),1)));

    Dpoints_gav = round(eeg_lat2point([0.3 1], [1 1], EEG.srate, [EEG.xmin EEG.xmax])); % index to time window of interest
    max_amp = max(max(GAV(:,Dpoints_gav(1):Dpoints_gav(2)))); % find the max amplitude within the time window of interest
    [chan lat_max_amp] = find(GAV == max_amp); % find the indicies for that amplitude
    lat_max_amp_ms = EEG.times(lat_max_amp); % convert the latency from sample point to ms / use this value while extracting ERPs
    GAVpeaks{g,2} = lat_max_amp_ms;
    GAVpeaks{g,3} = chan;

    % find midline peak electrode
    time_start_elecpeak = 400; % time window of interest (start)
    time_end_elecpeak = 600; % time winow of interest (end)
    trange_elecpeak = [time_start_elecpeak time_end_elecpeak]; % define the range
    Dpoints_elecpeak = round(eeg_lat2point(trange_elecpeak/1000, [1 1], EEG.srate, [EEG.xmin EEG.xmax])); % index to time window of interest

    % the list of channels is different when we use the mastoids as a reference, as opposed to the average
    % this is reflected in the number of channels; to avoid issues when
    % identifying the target electrode for the GAV peak, I
    % import the channel list of the corresponding reference-path here.
    % channels for avg ref is the same as for csd ref
    if g < 3
        chanlabels = crm; % channels reference mastoids
    else
        chanlabels = cra; % channels reference average (and csd)
    end

    for mp = 1:length(midline)
        mcis(mp) = find(strcmpi(midline{mp},chanlabels));
    end
    max_amp = max(max(GAV(mcis,Dpoints_elecpeak(1):Dpoints_elecpeak(2)))); % find the max amplitude within the time window of interest
    [mchan lat_max_amp] = find(GAV == max_amp);
    GAVpeaks{g,4} = chanlabels{mchan};
end

%% save peaks
% will be used in the next script
cd(storepath)
save('SAVpeaks.mat','SAVpeaks')
save('GAVpeaks.mat','GAVpeaks')