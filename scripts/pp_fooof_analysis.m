% pp_fooof_analysis.m
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
INPATH = fullfile(MAINPATH,'data\raw_data\pp_data_main_raw'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions');
addpath(FUNPATH);

%variables to edit
EPO_FROM = -0.8;
EPO_TILL = 0.7;
CHAN = 'Cz';
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};

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
        %load raw dataset (C_0001)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0001.vhdr'], [], []); %CHECK %channel locations already there
        %select channel
        EEG = pop_select( EEG, 'channel',{CHAN});
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
        %store dataset C_0001
        EEG.setname = [SUBJ '_talk_raw'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %load raw dataset (C_0005)
        EEG = pop_loadbv(fullfile(INPATH, SUBJ), ['av_' SUBJ '_C_0005.vhdr'], [], []); %CHECK %channel locations already there
        %select channel
        EEG = pop_select( EEG, 'channel',{CHAN});
        %rename events
        clear event
        for event = 1:length(EEG.event)
            if strcmp(EEG.event(event).type, 'S 64')
                EEG.event(event).type = 'listen';
            else
                continue
            end
        end
        %store dataset C_0005
        EEG.setname = [SUBJ '_listen_raw'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %merge datasets
        EEG = pop_mergeset( ALLEEG, [1 2], 0);   
        %epoching
        EEG = pop_epoch( EEG, {  }, [EPO_FROM  EPO_TILL], 'newname', 'Merged datasets epochs', 'epochinfo', 'yes');
        EEG.setname = [SUBJ '_talk_listen_raw'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %convert full data to fieldtrip format
        data = eeglab2fieldtrip(EEG, 'raw');

        %fooof analysis
        % Set up frequency analysis configuration with frequency range 20–60 Hz
        cfg               = [];
        cfg.foilim        = [20 60];    % analyze frequencies from 20 to 60 Hz
        cfg.pad           = 4;
        cfg.tapsmofrq     = 2;
        cfg.method        = 'mtmfft';
        cfg.output        = 'fooof_aperiodic';
        fractal = ft_freqanalysis(cfg, data);

        cfg.output        = 'pow';
        original = ft_freqanalysis(cfg, data);

        % Subtract the fractal component from the power spectrum
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2 - x1';
        oscillatory = ft_math(cfg, fractal, original);

        % Alternative definition of oscillatory component
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2 ./ x1';  % equivalent to 10^(log10(x2) - log10(x1))
        oscillatory_alt = ft_math(cfg, fractal, original);

        % Display the spectra on a linear frequency scale with limits set to [20 60]
        figure();
        subplot(1,2,1); hold on;
        plot(original.freq, log(original.powspctrm), 'k');
        plot(fractal.freq, log(fractal.powspctrm));
        plot(oscillatory.freq, log(oscillatory.powspctrm));
        xlim([20 60]);  % limit x-axis to frequencies between 20 and 60 Hz
        xlabel('Frequency (Hz)'); ylabel('log-power'); grid on;
        legend({'original', 'fractal', 'oscillatory = spectrum - fractal'}, 'location', 'southwest');
        title('mixed signal');

        subplot(1,2,2); hold on;
        plot(original.freq, log(original.powspctrm), 'k');
        plot(fractal.freq, log(fractal.powspctrm));
        plot(oscillatory_alt.freq, log(oscillatory_alt.powspctrm));
        xlim([20 60]);  % limit x-axis to frequencies between 20 and 60 Hz
        xlabel('Frequency (Hz)'); ylabel('log-power'); grid on;
        legend({'original', 'fractal', 'oscillatory = spectrum / fractal'}, 'location', 'southwest');
        title('oscillatory = spectrum / fractal');

        %sanity check variables
        subj_time = toc;
        ok_subj{subj,1} = SUBJ;
        ok_subj{subj,2} = subj_time;
    else
        marked_subj{end+1} = SUBJ;
    end
end

%display sanity check variables
marked_subj
ok_subj
check_done = 'OK'