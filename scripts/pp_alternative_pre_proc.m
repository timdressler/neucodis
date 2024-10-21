% pp_alternative_pre_proc.m
%
% Preproccessing pipeline for eeg-data, preceding a coherence analysis.
%
% Adapted from Moon et al. (2023)
% Tim Dressler, 11.09.2024

clear
close all
clc

%setup paths
SCRIPTPATH = cd;
%check if correct path is openend
if regexp(SCRIPTPATH, regexptranslate('wildcard','*neucodis\scripts')) == 1
    disp('Path OK')
else
    error('Path not OK')
end
MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH,"data\raw_data\pp_main_data_raw"); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_main_data_proc\pp_main_data_after_preproc_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
EVENTS = {'talk', 'listen'};
EPO_FROM = -0.8;
EPO_TILL = 0.7;
LCF = 1;
HCF = 60;
LCF_2 = 20;
HCF_2 = 55;
LCF_ICA = 1;
HCF_ICA = 30;
BL_FROM = -800;
BL_TILL= -600;
TF_FROM = -800;
TF_TILL = 600;
TF_BL_FROM = -600;
TF_BL_TILL = -400;
THRESH = 75;
SD_PROB = 3;
RESAM_ICA = 250;
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize marked subjects variable
MARKED_SUBJ = {};

%% load and merge data
for subj = 1:length(dircont_subj)
    %get current ID
    SUBJ = dircont_subj(subj).name;
    %check number of condition files
    dircont_cond1 = dir(fullfile(INPATH, [SUBJ '/*C_0001*.vhdr']));
    dircont_cond2 = dir(fullfile(INPATH, [SUBJ '/*C_0005*.vhdr']));
    if length(dircont_cond1) == 1 && length(dircont_cond2) == 1
        %start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %preprocessing dataset C_0001
        %load raw dataset (C_0001)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0001.vhdr'], [], []); %CHECK %channel locations already there
        %remove not needed channels
        EEG = pop_select( EEG, 'rmchannel',{'heogl','heogr','ml','mr', 'Lip'});
        %rename events
        for s = 2:length(EEG.event)
            if s < length(EEG.event)
                if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64') && ~strcmp(EEG.event(s+1).type, 'S 64')
                    EEG.event(s).type = 'talk';
                else
                    continue
                end
            elseif s == length(EEG.event)
                if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64')
                    EEG.event(s).type = 'talk';
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
        for s = 1:length(EEG.event)
            if strcmp(EEG.event(s).type, 'S 64')
                EEG.event(s).type = 'listen';
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
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
        EEG.setname = [SUBJ '_ICA_after_ICA'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

        %attach ICA weight to main data
        EEG = ALLEEG(3);
        CURRENTSET = 3;
        EEG = pop_editset(EEG,'run', [], 'icaweights','ALLEEG(8).icaweights', 'icasphere','ALLEEG(8).icasphere');
        %label ICA components with IC Label Plugin (Pion-Tonachini et al., 2019)
        EEG = pop_iclabel(EEG, 'default');
        EEG = pop_icflag(EEG, [0 0.2;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
        EEG = pop_subcomp( EEG, [], 0);
        %compute surface laplacian
        EEG.data=laplacian_perrinX(EEG.data,EEG.chanlocs.X,EEG.chanlocs.Y,EEG.chanlocs.Z);
        %20-55Hz bandpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',LCF_2,'hicutoff',HCF_2,'plotfreqz',0);
        %epoching
        EEG = pop_epoch( EEG, EVENTS, [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');
        %baseline removal
        EEG = pop_rmbase( EEG, [BL_FROM BL_TILL] ,[]);
        %threshold removal
        EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,-THRESH,THRESH,-0.8,0.699,2,1);
        EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,-THRESH,THRESH,-0.8,0.699,2,0);
        %probability-based removal
        EEG = pop_jointprob(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,SD_PROB,0,0,0,[],0);
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
        %end of preprocessing

        %save dataset
        EEG = pop_saveset(EEG, 'filename',[SUBJ '_after_preproc_proc.set'],'filepath', OUTPATH);
    else
        MARKED_SUBJ{end+1} = SUBJ;
    end
end

eeglab redraw






