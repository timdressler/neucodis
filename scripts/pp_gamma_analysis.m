% pp_gamma_analysis.m
%
% Plot ERSPs and Topoplots.
% Using preprocessed data.
% Preprocessing done with pp_coherence_pre_proc.m.
% Using talk & listen conditions.
% Plots ERSPs for each subject.
% Plots grand average ERSP & grand average Topoplot over all subjects.
%
% Tim Dressler, 02.11.2024

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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc_clean\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\gamma_analysis\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

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

%setup progress bar
wb = waitbar(0,'starting pp_gamma_analysis.m');

clear subj
for subj = 1:length(dircont_subj) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_coherence_preprocessed.set');
    %update progress bar
    waitbar(subj/length(dircont_subj),wb, [SUBJ ' pp_gamma_analysis.m'])
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
    xlim([-400 300])
    ylim([15 60])
    axis xy;
    colorbar;
    cb=colorbar; %CHECK %correct?
    cb.Title.String = "ERSP (dB)";
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title(['ERSP for Subject ' SUBJ]);

    set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
    saveas(gcf,fullfile(OUTPATH, ['subj_' SUBJ '_ersp.png']))

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

%% Grand Average Plots
%plot grand average ERSP
figure;
imagesc(ersp_times, ersp_freqs, ersp_GRANDAVERAGE);
axis xy;
xlim([-400 300])
ylim([15 60])
colorbar;
cb=colorbar; 
cb.Title.String = "ERSP (dB)";
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Grand Average ERSP');

set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'grand_average_ersp.png'))

%topoplot ERSP
figure;
tftopo(allersp_GRANDAVERAGE,alltimes(:,:,1),allfreqs(:,:,1), ...
    'timefreqs', [58 36; 60 38; -130 35; 60 43], 'chanlocs', EEG.chanlocs, 'showchan', chani)
sgtitle('Grand Average Topoplots')

set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'grand_average_topoplots.png'))

%% End of processing

%display sanity check variables
ok_subj = cell2table(ok_subj, 'VariableNames',{'subj','time'})
writetable(ok_subj,fullfile(OUTPATH, '_gamma_analysis_ok_subj.xlsx'))

check_done = 'OK'
save(fullfile(OUTPATH, '_gamma_analysis_data.mat'), 'check_done', 'allersp_GRANDAVERAGE', 'alltimes', 'allfreqs', 'EEG', 'chani')
save(fullfile(OUTPATH, '_gamma_analysis_plot_data.mat'), 'allersp_GRANDAVERAGE', 'alltimes', 'allfreqs', 'EEG', 'chani', ...
    'ersp_times', 'ersp_freqs', 'ersp_GRANDAVERAGE')


close(wb)
