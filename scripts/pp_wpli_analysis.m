% pp_wpli_analysis.m
%
% Description
%
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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_main_data_proc\pp_main_data_PREPROCESSED\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_main_data_proc\pp_main_analysis data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%selection of electrodes (left)
FRONTAL_L = {'F3', 'F5', 'F7', 'Fp1'};
TEMPORAL_L = {'T7', 'TP7', 'P7', 'FC5'};
OCCIPITAL_L = {'PO7', 'POz', 'Oz', 'O1'};

%check if number of electrodes is the same for each lobe (left)
switch length(FRONTAL_L) == length(TEMPORAL_L) && length(FRONTAL_L) == length(OCCIPITAL_L)
    case true
        disp('Electrodes OK')
        NUM_ELE = length(FRONTAL_L);
    otherwise
        error('Electrodes not OK')
end

%initialize electrode pairs variable (left)
PAIRS_L = {};

%generate electrode pairings (left)
r = 1;
for s = 1:NUM_ELE
    for v = 1:NUM_ELE
        PAIRS_L{r, 1} = FRONTAL_L{s};
        PAIRS_L{r, 2} = TEMPORAL_L{v};
        r = r+1;
    end
end

%initialize electrode pairs variable (control, left)
PAIRS_L_CONTROL = {};

%generate electrode pairings (control, left)
r = 1;
for s = 1:NUM_ELE
    for v = 1:NUM_ELE
        PAIRS_L_CONTROL{r, 1} = FRONTAL_L{s};
        PAIRS_L_CONTROL{r, 2} = OCCIPITAL_L{v};
        r = r+1;
    end
end

for subj = 1:length(dircont_subj)
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_PREPROCESSED.set');
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
    %check if ID matches dataset
    SUBJ_CHECK = strcmp(SUBJ, erase(EEG.setname, '_PREPROCESSED'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    %convert full data to fieldtrip format
    data = eeglab2fieldtrip(EEG, 'raw');
    %create datasets for each condition
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'listen'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [SUBJ '_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    data_listen = eeglab2fieldtrip(EEG, 'raw');
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'talk'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [SUBJ '_talk'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    data_talk = eeglab2fieldtrip(EEG, 'raw');

    %main analysis
    for elec_pair = 1:length(PAIRS_L) %loop over electrode pairs
        clear freq_talk freq_listen wpli_talk wpli_listen
        %TF transform via wavelets
        %setup
        cfg_freq = [];
        cfg_freq.method = 'wavelet';
        cfg_freq.output = 'powandcsd';
        cfg_freq.channel = {PAIRS_L{elec_pair,1}, PAIRS_L{elec_pair,2}};
        cfg_freq.keeptrials = 'yes';
        cfg_freq.foi = 25:1:50;
        cfg_freq.toi = -0.4:0.01:0.2;
        cfg_freq.width = 6;

        %TF analysis for talk condition
        freq_talk = ft_freqanalysis(cfg_freq, data_talk);
        %TF analysis for talk condition
        freq_listen = ft_freqanalysis(cfg_freq, data_listen);

        %connectivity analysis via debiased wPLI
        %setup
        cfg_conn = [];
        cfg_conn.method = 'wpli_debiased';
        cfg_conn.channelcmb = PAIRS_L;  % Paar 2: C3 - C4
        cfg_conn.keeptrials = 'yes';

        %connectivity analysis for talk condition
        wpli_talk = ft_connectivityanalysis(cfg_conn, freq_talk);
        wpli_talk_extracted = squeeze(wpli_talk.wpli_debiasedspctrm);
        %connectivity analysis for listen condition
        wpli_listen = ft_connectivityanalysis(cfg_conn, freq_listen);
        wpli_listen_extracted = squeeze(wpli_listen.wpli_debiasedspctrm);

        %store wPLI values over all pairs
        wpli_talk_extracted_ALL_PAIRS(:,:,elec_pair) = wpli_talk_extracted;
        wpli_listen_extracted_ALL_PAIRS(:,:,elec_pair) = wpli_listen_extracted;
    end

    %get mean wPLI over all pairs
    wpli_talk_extracted_AVERAGE_PAIRS = mean(wpli_talk_extracted_ALL_PAIRS,3);
    wpli_listen_extracted_AVERAGE_PAIRS = mean(wpli_listen_extracted_ALL_PAIRS,3);
    %store the average values in original structure
    wpli_talk_AVERAGE = wpli_listen;
    wpli_talk_AVERAGE.wpli_debiasedspctrm = wpli_talk_extracted_AVERAGE_PAIRS;
    wpli_talk_AVERAGE.label = {'FRONT_TEMP'};
    wpli_listen_AVERAGE = wpli_listen;
    wpli_listen_AVERAGE.wpli_debiasedspctrm = wpli_listen_extracted_AVERAGE_PAIRS;
    wpli_listen_AVERAGE.label = {'FRONT_TEMP'};

    %store structures in cell
    wpli_talk_AVERAGE_ALL_SUBJ{subj} = wpli_talk_AVERAGE;
    wpli_listen_AVERAGE_ALL_SUBJ{subj} = wpli_listen_AVERAGE;

    %sanity checks
    subj_time = toc;
    OK_SUBJ{subj,1} = SUBJ;
    OK_SUBJ{subj,2} = SUBJ_CHECK;
    OK_SUBJ{subj,3} = subj_time;
end

OK_SUBJ

%get grand mean over all subjects
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = 'wpli_debiasedspctrm';      
wpli_talk_GRANDAVERAGE = ft_freqgrandaverage(cfg, wpli_talk_AVERAGE_ALL_SUBJ{:});
wpli_listen_GRANDAVERAGE = ft_freqgrandaverage(cfg, wpli_listen_AVERAGE_ALL_SUBJ{:});




