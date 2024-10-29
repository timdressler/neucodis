% pp_erp_pre_proc.m
%
% Preproccessing pipeline for eeg-data, preceding an erp analysis.
%
% Adapted from Moon et al. (2023)
% Tim Dressler, 11.09.2024

clear
close all
clc

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
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions');
addpath(FUNPATH);

%variables to edit
EPO_FROM = -0.25;
EPO_TILL = 0.750;
EPO_FROM_BL = -2;
EPO_TILL_BL = 1;
LCF = 1;
HCF = 60;
LCF_2 = 1;
HCF_2 = 30;
LCF_ICA = 1;
HCF_ICA = 30;
BL_FROM = -250;
THRESH = 75;
SD_PROB = 3;
RESAM_ICA = 250;
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};
all_epo_lats_OK_ALL = [];

%% load and merge data
clear subj
for subj = 1:length(dircont_subj)
    %get current ID
    SUBJ = dircont_subj(subj).name;
    %check number of condition files
    dircont_cond1 = dir(fullfile(INPATH, [SUBJ '/*C_0001*.vhdr']));
    dircont_cond2 = dir(fullfile(INPATH, [SUBJ '/*C_0005*.vhdr']));
    if length(dircont_cond1) == 1 && length(dircont_cond2) == 1
        tic;
        %start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %preprocessing dataset C_0001
        %load raw dataset (C_0001)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0001.vhdr'], [], []); %CHECK %channel locations already there
        %remove not needed channels
        EEG = pop_select( EEG, 'rmchannel',{'heogl','heogr','ml','mr', 'Lip'});
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
        %1-60Hz bandpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',HCF,'plotfreqz',0);
        %run PREP pipeline (Bigdely-Shamlo et al., 2015)
        %settings
        params = struct();
        params.lineFrequencies = 50:50:EEG.srate/2-50; %set line noise to 50Hz
        params.ignoreBoundaryEvents = true;  %ingore boundary events
        params.reportMode = 'skipReport';  %suppress report
        params.keepFiltered = true; %remove trend
        %run
        EEG = pop_prepPipeline(EEG, params);
        %store dataset C_0001
        EEG.setname = [SUBJ '_talk_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %preprocessing dataset C_0005
        %load raw dataset (C_0005)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0005.vhdr'], [], []); %CHECK %channel locations already there
        %remove not needed channels
        EEG = pop_select( EEG, 'rmchannel',{'heogl','heogr','ml','mr', 'Lip'});
        %rename events
        clear event
        for event = 1:length(EEG.event)
            if strcmp(EEG.event(event).type, 'S 64')
                EEG.event(event).type = 'listen';
            else
                continue
            end
        end
        %1-60Hz bandpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',HCF,'plotfreqz',0);
        %run PREP pipeline (Bigdely-Shamlo et al., 2015)
        %settings
        params = struct();
        params.lineFrequencies = 50:50:EEG.srate/2-50; %set line noise to 50Hz
        params.ignoreBoundaryEvents = true;  %ingore boundary events
        params.reportMode = 'skipReport';  %suppress report
        params.keepFiltered = true; %remove trend
        %run
        EEG = pop_prepPipeline(EEG, params);
        %store dataset C_0005
        EEG.setname = [SUBJ '_listen_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %preprocessing on merged datasets
        %merge datasets
        EEG = pop_mergeset( ALLEEG, [1 2], 0);
        EEG.setname = [SUBJ '_merged_after_PREP'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %ICA preprocessing
        %reload and merge raw data
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0001.vhdr'], [], []); %CHECK %channel locations already there
        EEG.setname = [SUBJ '_talk_raw_before_ICA'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0005.vhdr'], [], []); %CHECK %channel locations already there
        EEG.setname = [SUBJ '_listen_raw_before_ICA'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG = pop_mergeset( ALLEEG, [1 2], 0);
        EEG.setname = [SUBJ '_merged_raw_before_ICA'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %1-30Hz bandpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',LCF_ICA,'hicutoff',HCF_ICA,'plotfreqz',0);
        %resample to 250Hz
        EEG = pop_resample( EEG, RESAM_ICA);
        EEG.setname = [SUBJ '_ICA_after_preproc'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %run ICA
        % % EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
        EEG.setname = [SUBJ '_ICA_after_ICA'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %attach ICA weight to main data
        EEG = ALLEEG(3);
        CURRENTSET = 3;
        % % EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(8).icaweights', 'icasphere','ALLEEG(8).icasphere');
        %label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
        % % EEG = pop_iclabel(EEG, 'default');
        % % EEG = pop_icflag(EEG, [0 0.2;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
        % % EEG = pop_subcomp( EEG, [], 0);
        %20-55Hz bandpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',LCF_2,'hicutoff',HCF_2,'plotfreqz',0);
        %rereference to average of mastoids
        EEG = pop_reref( EEG, [31 32] ); %CHECK %schon average?
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
                error('ERROR: problem with epoch latency variable, more that 2 or less than one event per epoch')
        end

        %baseline removal 
        %remove baseline from BL_FROM until trial onset
        %trial onset for talk = fixation cross
        %trial onset for listen = playback onset
        clear epo
        for epo = 1:length(EEG.epoch)
            %get latencies
            trial_onset = EEG.epoch(epo).eventlatency{1};
            bl_start_lat = trial_onset+BL_FROM;
            bl_end_lat = trial_onset;
            %convert latencies to samples
            [~, bl_start_sam] = min(abs(EEG.times-bl_start_lat));
            [~, bl_end_sam] = min(abs(EEG.times-bl_end_lat));
            %get and subtract baseline
            bl = mean(EEG.data(:,bl_start_sam:bl_end_sam,epo), 2);
            EEG.data(:,:,epo) = EEG.data(:,:,epo) - bl;
        end
        %epoching
        EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');
        %threshold removal
        EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,-THRESH,THRESH,EPO_FROM,EPO_TILL,0,0);
        %probability-based removal
        EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
        %end of preprocessing

        %save dataset
        EEG.setname = [SUBJ '_erp_preprocessed'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG = pop_saveset(EEG, 'filename',[SUBJ '_erp_preprocessed.set'],'filepath', OUTPATH);
        
        %sanity check variables
        subj_time = toc;
        ok_subj{subj,1} = SUBJ;
        ok_subj{subj,2} = subj_time;
    else
        marked_subj{end+1} = SUBJ;
    end
    all_epo_lats_OK_ALL(end+1) = all_epo_lats_OK;
end

%display sanity check variables
marked_subj
ok_subj
all(all_epo_lats_OK_ALL)
check_done = 'OK'

%%eeglab redraw