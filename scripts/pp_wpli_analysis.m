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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_main_data_proc\pp_main_data_after_preproc_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_main_data_proc\pp_main_analysis data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%selection of electrodes (left)
FRONTAL_L = {'F3', 'F5', 'F7', 'F9'};
TEMPORAL_L = {'T7', 'TP9', 'T3', 'FC5'};
OCCIPITAL_L = {'O9', 'O1', 'PO9', 'PO7'};

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
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
    %convert data to fieldtrip format
    data = eeglab2fieldtrip(EEG, 'raw');
    %create datasets for each condition
    EEG = pop_selectevent( EEG, 'latency','-2<=2','type',{'listen'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    data_talk = eeglab2fieldtrip(EEG, 'raw');
    EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',{'talk'},'deleteevents','off','deleteepochs','on','invertepochs','off');
    data_listen = eeglab2fieldtrip(EEG, 'raw');

    %main analysis
    %TF transform via wavelets
    %setup
    cfg_freq = [];
    cfg_freq.method = 'wavelet';
    cfg_freq.output = 'powandcsd';
    cfg_freq.channel = {'F7', 'T7'};
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
    cfg_conn.keeptrials = 'yes';

    %connectivity analysis for talk condition
    wpli_talk = ft_connectivityanalysis(cfg_conn, freq_talk);
    %connectivity analysis for listen condition
    wpli_listen = ft_connectivityanalysis(cfg_conn, freq_listen);






end

eeglab redraw


