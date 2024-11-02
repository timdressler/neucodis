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
%selection of electrodes (left)
FRONTAL_L = {'F3', 'F5', 'F7', 'Fp1'};
TEMPORAL_L = {'T7', 'TP7', 'P7', 'FC5'};
OCCIPITAL_L = {'PO7', 'POz', 'Oz', 'O1'};
%selection of electrodes (right)
FRONTAL_R = {'F2', 'F6', 'F8', 'Fp2'};
TEMPORAL_R = {'T8', 'TP8', 'P8', 'FC6'};
OCCIPITAL_R = {'PO8', 'POz', 'Oz', 'O2'};

%sanity check
%check if number of electrodes is the same for each lobe (left)
switch length(FRONTAL_L) == length(TEMPORAL_L) && length(FRONTAL_L) == length(OCCIPITAL_L)
    case true
        disp('Electrodes OK')
        num_ele = length(FRONTAL_L);
    otherwise
        error('Electrodes not OK')
end

%setup electrode pairs
all_pairs(1).name = 'Fronto-Temporal_L';
all_pairs(2).name = 'Fronto-Temporal_R';
all_pairs(3).name = 'Fronto-Occipital_L';
all_pairs(4).name = 'Fronto-Occipital_R';

all_pairs(1).elec1 = FRONTAL_L;
all_pairs(1).elec2 = TEMPORAL_L;
all_pairs(2).elec1 = FRONTAL_R;
all_pairs(2).elec2 = TEMPORAL_R;
all_pairs(3).elec1 = FRONTAL_L;
all_pairs(3).elec2 = OCCIPITAL_L;
all_pairs(4).elec1 = FRONTAL_R;
all_pairs(4).elec2 = OCCIPITAL_R;

for pairs = 1:length(all_pairs)
    r = 1;
    for elec1 = 1:num_ele
        for elec2 = 1:num_ele
            all_pairs(pairs).pairs{r, 1} = all_pairs(pairs).elec1{elec1};
            all_pairs(pairs).pairs{r, 2} = all_pairs(pairs).elec2{elec2};
            r = r+1;
        end
    end
end

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};

clear pairs
r = 1;
for pairs = 1:length(all_pairs) %loop over electode pair-sets
    r_start = r;
    clear subj subj_time subj_check
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
        %rename dataset
        EEG.setname = [SUBJ '_talk_listen'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        
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
        clear elec_pair
        for elec_pair = 1:length(all_pairs(pairs).pairs) %loop over electrode pairs
            clear freq_talk freq_listen wpli_talk wpli_listen
            %TF transform via wavelets
            %setup
            cfg_freq = [];
            cfg_freq.method = 'wavelet';
            cfg_freq.output = 'powandcsd';
            cfg_freq.channel = {all_pairs(pairs).pairs{elec_pair,1}, all_pairs(pairs).pairs{elec_pair,2}};
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
            % % cfg_conn.channelcmb = PAIRS_L;  % Paar 2: C3 - C4
            cfg_conn.keeptrials = 'yes';

            %connectivity analysis for talk condition
            wpli_talk = ft_connectivityanalysis(cfg_conn, freq_talk);
            %baseline correction

            wpli_talk_extracted = squeeze(wpli_talk.wpli_debiasedspctrm);
            %connectivity analysis for listen condition
            wpli_listen = ft_connectivityanalysis(cfg_conn, freq_listen);
            %baseline correction
            
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
        wpli_talk_AVERAGE.wpli_debiasedspctrm(1,:,:) = wpli_talk_extracted_AVERAGE_PAIRS;
        wpli_talk_AVERAGE.label = {all_pairs(pairs).name};
        wpli_talk_AVERAGE.dimord = 'chan_freq_time';
        wpli_listen_AVERAGE = wpli_listen;
        wpli_listen_AVERAGE.wpli_debiasedspctrm(1,:,:) = wpli_listen_extracted_AVERAGE_PAIRS;
        wpli_listen_AVERAGE.label = {all_pairs(pairs).name};
        wpli_listen_AVERAGE.dimord = 'chan_freq_time';


        %store structures in cell
        wpli_talk_AVERAGE_ALL_SUBJ{subj} = wpli_talk_AVERAGE;
        wpli_listen_AVERAGE_ALL_SUBJ{subj} = wpli_listen_AVERAGE;

        %sanity checks
        subj_time = toc;
        ok_subj{r,1} = SUBJ;
        ok_subj{r,2} = subj_check;
        ok_subj{r,3} = subj_time;
        r = r+1;
        if r <= length(dircont_subj)*4;
            r_end = r;
        end
    end


    %get grand average over all subjects
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'wpli_debiasedspctrm';
    wpli_talk_GRANDAVERAGE = ft_freqgrandaverage(cfg, wpli_talk_AVERAGE_ALL_SUBJ{:});
    wpli_listen_GRANDAVERAGE = ft_freqgrandaverage(cfg, wpli_listen_AVERAGE_ALL_SUBJ{:});

    %store cells in structure
    all_wpli(pairs).name = {all_pairs(pairs).name};
    all_wpli(pairs).talk_GA = wpli_talk_GRANDAVERAGE;
    all_wpli(pairs).listen_GA = wpli_listen_GRANDAVERAGE;

    % % cfg = [];
    % % cfg.method           = 'montecarlo';
    % % cfg.statistic        = 'ft_statfun_depsamplesT';
    % % cfg.correctm         = 'cluster';
    % % cfg.clusteralpha     = 0.05;
    % % cfg.clusterstatistic = 'maxsum';
    % % cfg.alpha            = 0.025;
    % % cfg.numrandomization = 500;

    % % all_wpli(pairs).comparison = ft_freqstatistics(cfg,all_wpli(pairs).talk_GA,all_wpli(pairs).listen_GA);


    %sanity checks
    for temp = r_start:r_end
        ok_subj{temp,4} = all_wpli(pairs).name{1};
    end
end

%display sanity check variables
ok_subj
check_done = 'OK'

%% Plot wPLI
time_vector = wpli_listen.time;
freq_vector = wpli_listen.freq;
figure;
imagesc(time_vector, freq_vector, squeeze(mean(wpli_listen.wpli_debiasedspctrm, 1)));
axis xy;
xlabel('Zeit (s)');
ylabel('Frequenz (Hz)');
title('wPLI Heatmap über Zeit und Frequenz (Wavelet) (listen)');
colorbar;
caxis([0, max(max(squeeze(mean(wpli_listen.wpli_debiasedspctrm, 1))))]);

time_vector = wpli_talk.time;
freq_vector = wpli_talk.freq;
figure;
imagesc(time_vector, freq_vector, squeeze(mean(wpli_talk.wpli_debiasedspctrm, 1)));
axis xy;
xlabel('Zeit (s)');
ylabel('Frequenz (Hz)');
title('wPLI Heatmap über Zeit und Frequenz (Wavelet) (talk)');
colorbar;
caxis([0, max(max(squeeze(mean(wpli_talk.wpli_debiasedspctrm, 1))))]);

