% pp_wpli_analysis.m
%
% Description
%
% Tim Dressler, 11.09.2024

clear
close all
clc

%critical variables to edit
%select whether to use simulated or real data
SIM_DATA = 1; %1 for simulated data %0 for real data

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
switch SIM_DATA
    case 0
        INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
        disp('Real data used!')
    case 1
        INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_simulated\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
        warning('Simulated data used!')
    otherwise
        error('error')
end
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
WPLI_BL_FROM = -0.4;
WPLI_BL_TILL = -0.3;
TF_F_FROM = 25;
TF_F_TILL = 50;
TF_T_FROM = -0.4;
TF_T_TILL = 0.2;
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
for pairs = 1:length(all_pairs) %loop over electrode pair-sets
    r_start = r;
    clear subj subj_time subj_check
    for subj = 1:length(dircont_subj) %loop over subjects
        tic;
        %get current ID
        switch SIM_DATA
            case 0
                SUBJ = erase(dircont_subj(subj).name, '_coherence_preprocessed.set');
            case 1
                SUBJ = erase(dircont_subj(subj).name, '.set');
        end

        %import data
        %start eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %import dataset
        EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
        %sanity check
        %check if ID matches dataset
        switch SIM_DATA
            case 0
                subj_check = strcmp(SUBJ, erase(EEG.setname, '_coherence_preprocessed'));
            case 1
                subj_check = strcmp(SUBJ, EEG.setname);
        end
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
            cfg_freq.foi = TF_F_FROM:1:TF_F_TILL;
            cfg_freq.toi = TF_T_FROM:0.01:TF_T_TILL;
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
            %extract wpli values
            wpli_talk_extracted = squeeze(wpli_talk.wpli_debiasedspctrm);
            %convert latencies to samples
            [~, wpli_bl_start_sam] = min(abs(wpli_talk.time-WPLI_BL_FROM));
            [~, wpli_bl_end_sam] = min(abs(wpli_talk.time-WPLI_BL_TILL));
            %baseline correction
            % % wpli_talk_extracted = wpli_talk_extracted - mean(wpli_talk_extracted(:,wpli_bl_start_sam:wpli_bl_end_sam),2);
            %connectivity analysis for listen condition
            wpli_listen = ft_connectivityanalysis(cfg_conn, freq_listen);
            %extract wpli values
            wpli_listen_extracted = squeeze(wpli_listen.wpli_debiasedspctrm);
            %convert latencies to samples
            [~, wpli_bl_start_sam] = min(abs(wpli_talk.time-WPLI_BL_FROM));
            [~, wpli_bl_end_sam] = min(abs(wpli_talk.time-WPLI_BL_TILL));
            %baseline correction
            % % wpli_talk_extracted = wpli_talk_extracted - mean(wpli_talk_extracted(:,wpli_bl_start_sam:wpli_bl_end_sam),2);
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

    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic         = 'depsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
    cfg.correctm          = 'cluster';
    cfg.clusteralpha      = 0.05;                   % alpha level of the sample-specific test statistic that will be used for thresholding
    cfg.clustertail       = 0;
    cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
    cfg.alpha             = 0.05;                   % alpha level of the permutation test
    cfg.numrandomization  = 1000;                   % number of draws from the permutation distribution
    cfg.ivar              = 1;                      % the index of the independent variable in the design matrix
    cfg.neighbours        = [];                     % there are no spatial neighbours, only in time and frequency
    cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
    cfg.parameter = 'wpli_debiasedspctrm';

    subj_design = length(dircont_subj);
    design = zeros(2,2*subj_design);
    for i = 1:subj_design
        design(1,i) = i;
    end
    for i = 1:subj_design
        design(1,subj_design+i) = i;
    end
    design(2,1:subj_design)        = 1;
    design(2,subj_design+1:2*subj_design) = 2;

    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sanity check
    %check dimensions of wPLI structure
    clear ok_dim
    switch all(size(all_wpli(pairs).talk_GA.wpli_debiasedspctrm) == size(all_wpli(pairs).listen_GA.wpli_debiasedspctrm),'all')
        case true
            ok_dim = 1;
        otherwise
            ok_dim = 0;
    end

    %cluster test
    all_wpli(pairs).comparison = ft_freqstatistics(cfg,all_wpli(pairs).talk_GA,all_wpli(pairs).listen_GA);

    %sanity checks
    clear temp
    for temp = r_start:r_end
        ok_subj{temp,4} = all_wpli(pairs).name{1};
        ok_subj{temp,5} = ok_dim;
    end
end

%display sanity check variables
ok_subj
check_done = 'OK'
switch SIM_DATA
    case 1
        warning('Simulated data has been used')
    case 0
        disp('Real data has been used')
end

%% Plot wPLI
set(0,'DefaultTextInterpreter','none') 
close all

clear pairs
for pairs = 1:length(all_pairs)
    %plot wPLI listen condition
    time_vector = all_wpli(pairs).listen_GA.time;
    freq_vector = all_wpli(pairs).listen_GA.freq;
    figure;
    imagesc(time_vector, freq_vector, squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1)));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['wPLI values for listen condition - ' all_wpli(pairs).name]);
    colorbar;
    caxis([0 0.5]);

    %plot wPLI talk condition
     time_vector = all_wpli(pairs).talk_GA.time;
    freq_vector = all_wpli(pairs).talk_GA.freq;
    figure;
    imagesc(time_vector, freq_vector, squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['wPLI values for talk condition - ' all_wpli(pairs).name]);
    colorbar;
    caxis([0 0.5]);

    %plot wPLI difference between talk and listen condition
     time_vector = all_wpli(pairs).listen_GA.time;
    freq_vector = all_wpli(pairs).listen_GA.freq;
    figure;
    imagesc(time_vector, freq_vector,squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1)) - ...
        squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['wPLI difference values for talk vs. listen condition - ' all_wpli(pairs).name]);
    colorbar;
    caxis([-0.3 0.3]);
end

%% Plot Cluster-Based Permutation Test


time = all_wpli(1).talk_GA.time;           % Time vector
freq = all_wpli(1).talk_GA.freq;           % Frequency vector
sig_clusters = squeeze(all_wpli(1).comparison.mask);  % Significance mask
talk_GA_extracted = squeeze(mean(all_wpli(1).talk_GA.wpli_debiasedspctrm,1));
listen_GA_extracted = squeeze(mean(all_wpli(1).listen_GA.wpli_debiasedspctrm,1));
effect = talk_GA_extracted - talk_GA_extracted;

% Plot the TF map
figure;
imagesc(time, freq, effect);
axis xy;  % Ensure the y-axis (frequency) is oriented correctly
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('wPLI Difference Time-Frequency Plot with Significant Clusters');

% Customize colormap and color limits for wPLI
caxis([-0.5, 0.5]);  % Adjust this based on wPLI difference range
colormap("parula");   % Using a perceptually uniform colormap

hold on;

% Overlay significant clusters using contour
contour(time, freq, sig_clusters, 1, 'LineColor', 'r', 'LineWidth', 2);

hold off;


%%%

any(all_wpli(1).comparison.mask, 'all')
disp(any(all_wpli(1).comparison.mask(:)))

