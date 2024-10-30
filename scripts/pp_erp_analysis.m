% pp_erp_analysis.m
%
% Description
%
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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};

for subj = 1:length(dircont_subj) %loop over subjects
        tic;
        %get current ID
        SUBJ = erase(dircont_subj(subj).name, '_erp_preprocessed.set');
        %import data
        %start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %import dataset
        EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
        %sanity check
        %check if ID matches dataset
        subj_check = strcmp(SUBJ, erase(EEG.setname, '_erp_preprocessed'));
        %get ERP
        ERP = mean(EEG.data, 3);


        %sanity checks
        subj_time = toc;
        ok_subj{subj,1} = SUBJ;
        ok_subj{subj,2} = subj_check;
        ok_subj{subj,3} = subj_time;
end




ok_subj
























