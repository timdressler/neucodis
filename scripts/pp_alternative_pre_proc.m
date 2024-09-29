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
MAINPATH = "C:/Users/timdr/OneDrive/Uni_Oldenburg/3_Semester\Module/Pratical_Project/Analysis";
INPATH = fullfile(MAINPATH,"data/raw_data/pp_main_data_raw/");
addpath("C:/Users/timdr/OneDrive/Uni_Oldenburg/3_Semester\Module/Pratical_Project/Analysis/neucodis/functions")

%% load and merge data
%start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%load raw dataset (C_0001)
EEG = pop_loadbv(fullfile(INPATH, 'P136'), 'av_P136_C_0001.vhdr', [], []);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%load raw dataset (C_0005)
EEG = pop_loadbv(fullfile(INPATH, 'P136'), 'av_P136_C_0005.vhdr', [], []);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%merge datasets
%%EEG = pop_mergeset( ALLEEG, [1  2], 0);

%% filter the data 
%1-60Hz bandpass
EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',60,'plotfreqz',0);





eeglab redraw