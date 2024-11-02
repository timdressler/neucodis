% pp_gamma_analysis.m
%
% Description
%
% Tim Dressler, 02.11.2024

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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
CHAN = 'Cz';
EVENTS = {'talk', 'listen'};
TF_FROM = -400;
TF_TILL = 500;
TF_BL_FROM = -200;
TF_BL_TILL = -100;

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
ok_subj = {};

for subj = 1:length(dircont_subj) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_coherence_preprocessed.set');
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
    %sanity check
    %check if ID matches dataset
    subj_check = strcmp(SUBJ, erase(EEG.setname, '_coherence_preprocessed'));
    %get channel ID
    chani = find(strcmp({EEG.chanlocs.labels}, CHAN));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    %compute ERSP
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = pop_newtimef(EEG, ...
        1, chani, [EEG.xmin EEG.xmax]*1000, [3 0.5], 'maxfreq', 60, 'padratio', 16, ...
        'plotphase', 'off', 'timesout', 60, 'alpha', .05, 'plotersp','off', 'plotitc','off');
    %store ERSP values
    all_ersp{subj} = ersp;
    %store time and frequency vectors
    ersp_times = times;
    ersp_freqs = freqs;

    %plot ERSP for each subject
    figure;
    imagesc(ersp_times, ersp_freqs, ersp);
    axis xy;
    colorbar;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title(['ERSP for Subject ' SUBJ]);

    %compute ERSPs over all electrodes TF-topoplots
    for elec = 1:EEG.nbchan
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = pop_newtimef(EEG, ...
            1, elec, [EEG.xmin EEG.xmax]*1000, [3 0.5], 'maxfreq', 60, 'padratio', 16, ...
            'plotphase', 'off', 'timesout', 60, 'alpha', .05, 'plotersp','off', 'plotitc','off');
        %create empty arrays if first electrode
        if elec == 1
            allersp = zeros([ size(ersp) EEG.nbchan]);
            allitc = zeros([ size(itc) EEG.nbchan]);
            allpowbase = zeros([ size(powbase) EEG.nbchan]);
            alltimes = zeros([ size(times) EEG.nbchan]);
            allfreqs = zeros([ size(freqs) EEG.nbchan]);
            allerspboot = zeros([ size(erspboot) EEG.nbchan]);
            allitcboot = zeros([ size(itcboot) EEG.nbchan]);
        end
        allersp (:,:,elec) = ersp;
        allitc (:,:,elec) = itc;
        allpowbase (:,:,elec) = powbase;
        alltimes (:,:,elec) = times;
        allfreqs (:,:,elec) = freqs;
        allerspboot (:,:,elec) = erspboot;
        allitcboot (:,:,elec) = itcboot;
    end

    all_allersp(:,:,:,subj) = allersp;

    %sanity check variables
    subj_time = toc;
    ok_subj{subj,1} = SUBJ;
    ok_subj{subj,2} = subj_time;
end

%re-arrange ERSP values
all_ersp = cat(3, all_ersp{:});

%get grand average over all subjects
ersp_GRANDAVERAGE = mean(all_ersp, 3);
allersp_GRANDAVERAGE = mean(all_allersp,4);

%plot grand average ERSP
figure;
imagesc(ersp_times, ersp_freqs, ersp_GRANDAVERAGE);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Grand Average ERSP');

%topoplot ERSP
figure;
tftopo(allersp_GRANDAVERAGE,alltimes(:,:,1),allfreqs(:,:,1), ...
    'timefreqs', [58 36; 70 48; 70 38; 60 43], 'chanlocs', EEG.chanlocs, 'showchan', chani)

%display sanity check variables
ok_subj
check_done = 'OK'
