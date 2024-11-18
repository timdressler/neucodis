clear
close all
clc

set(0,'DefaultTextInterpreter','none')

%setup paths
SCRIPTPATH = cd;

%sanity check 
%check if paths are correct
if regexp(SCRIPTPATH, regexptranslate('wildcard','*neucodis\scripts')) == 1
    disp('Path OK')
else
    error('Path not OK')
end

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH,'data\raw_data\pp_data_main_raw'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH)


%variables to edit
EPO_FROM_BL = -2;
EPO_TILL_BL = 1;
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};
all_epo_lats_OK_ALL = [];

%setup progress bar
wb = waitbar(0,'starting pp.m');

clear subj
for subj = 1:length(dircont_subj)
    %get current ID
    SUBJ = dircont_subj(subj).name;
    %update progress bar
    waitbar(subj/length(dircont_subj),wb, [SUBJ ' pp.m'])
    %check number of condition files
    dircont_cond1 = dir(fullfile(INPATH, [SUBJ '/*C_0001*.vhdr']));
    dircont_cond2 = dir(fullfile(INPATH, [SUBJ '/*C_0005*.vhdr']));
    if length(dircont_cond1) == 1 && length(dircont_cond2) == 1
        tic;
        %start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %preprocessing dataset C_0001
        %load raw dataset (C_0001)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0001.vhdr'], [], []); 
        %rename events
        clear event
        for event = 2:length(EEG.event)
            if event < length(EEG.event)
                if strcmp(EEG.event(event).type, 'S 64') && ~strcmp(EEG.event(event-1).type, 'S 64') && ~strcmp(EEG.event(event+1).type, 'S 64')
                    EEG.event(event).type = 'talk';
                else
                    continue
                end
            elseif event == length(EEG.event)
                if strcmp(EEG.event(event).type, 'S 64') && ~strcmp(EEG.event(event-1).type, 'S 64')
                    EEG.event(event).type = 'talk';
                else
                    continue
                end
            end
        end
        %store dataset C_0001
        EEG.setname = [SUBJ '_talk_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %preprocessing dataset C_0005
        %load raw dataset (C_0005)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0005.vhdr'], [], []); %CHECK %channel locations already there
        %rename events
        clear event
        for event = 1:length(EEG.event)
            if strcmp(EEG.event(event).type, 'S 64')
                EEG.event(event).type = 'listen';
            else
                continue
            end
        end
        %store dataset C_0005
        EEG.setname = [SUBJ '_listen_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %preprocessing on merged datasets
        %merge datasets
        EEG = pop_mergeset( ALLEEG, [1 2], 0);
        EEG.setname = [SUBJ '_merged_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %baseline epoching 
        EEG = pop_epoch( EEG, EVENTS, [EPO_FROM_BL        EPO_TILL_BL], 'epochinfo', 'yes');
    
        %sanity checks before baseline removal
        %check whether all events contain two (fixation cross + speech
        %onset) events (for talk) or one (playback onset) event (for
        %listen)
        %check whether the second (for talk) or the first (for listen)
        %event latency is zero
        clear epo all_epo_lats_OK epo_lats
        epo_lats = {EEG.epoch.eventlatency};
        all_epo_lats_OK = all(cellfun(@(epo) (isequal(size(epo), [1, 1]) && epo{1} == 0) || ...
            (isequal(size(epo), [1, 2]) && epo{2} == 0) && epo{1} < epo{2} && epo{1} > -2000, epo_lats));
        switch all_epo_lats_OK
            case true
                disp('epoch latency variable OK')
            case false
                warning([SUBJ ': problem with epoch latency variable, more that 2 or less than one event per epoch'])
        end

        %get latencies
        clear epo
        for epo = 1:length(EEG.epoch)
            %get latencies
            lats(epo) = EEG.epoch(epo).eventlatency{1}; 
        end
        max_lat(subj) = min(lats);
        %end of preprocessing

        %sanity check variables
        subj_time = toc;
        ok_subj{subj,1} = SUBJ;
        ok_subj{subj,2} = subj_time;
    else
        marked_subj{end+1} = SUBJ;
    end
    all_epo_lats_OK_ALL(end+1) = all_epo_lats_OK;
end

%% End of processing



close(wb)


