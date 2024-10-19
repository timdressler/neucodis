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




    

    %main analysis
%% Frequenzanalyse zur Berechnung der Kreuzspektraldichte
cfg_freq = [];
cfg_freq.method = 'mtmconvol';               % Zeit-Frequenz-Analyse verwenden
cfg_freq.output = 'powandcsd';               % Power und Kreuzspektraldichte berechnen
cfg_freq.channel = {'F7', 'T7'};             % Die gewünschten Kanäle
cfg_freq.keeptrials = 'yes';                 % Trials beibehalten
cfg_freq.taper = 'dpss';                      % DPSS-Taper verwenden
cfg_freq.tapsmofrq = 7.5;                       % Glättungsparameter anpassen
cfg_freq.foi = 1:1:50;                        % Frequenzen von 1 Hz bis 30 Hz, um die Datenlänge nicht zu überschreiten
cfg_freq.t_ftimwin = 5 ./ cfg_freq.foi;      % Zeitfenster für jede Frequenz
cfg_freq.toi = data.time{1};                  % Zeitpunkte basierend auf den Epochen
freq_data = ft_freqanalysis(cfg_freq, data);  % Frequenzanalyse durchführen

%% WPLI-Analyse
cfg_conn = [];
cfg_conn.method = 'wpli_debiased';            % Debiased wPLI verwenden
cfg_conn.keeptrials = 'yes';                  % Trials beibehalten
cfg_conn.channelcmb = {'F7', 'T7'};           % Kanal-Kombination
wpli = ft_connectivityanalysis(cfg_conn, freq_data); % WPLI berechnen

%% Plotten der wPLI-Werte
%% Heatmap der wPLI-Werte
figure;

% Zeit- und Frequenzvektoren
time_vector = wpli.time;  % Zeitpunkte für den Plot
freq_vector = wpli.freq;   % Frequenzen

% Heatmap der wPLI-Werte
imagesc(time_vector, freq_vector, squeeze(mean(wpli.wpli_debiasedspctrm, 1))); 
axis xy; % Achse umkehren
xlabel('Zeit (s)');
ylabel('Frequenz (Hz)');
title('wPLI Heatmap über Zeit und Frequenz');
colorbar; % Farbskala anzeigen
caxis([0, max(max(squeeze(mean(wpli.wpli_debiasedspctrm, 1))))]); % Farbskala anpassen




end

eeglab redraw



%% Frequenzanalyse zur Berechnung der Kreuzspektraldichte mit Wavelet-Transformation
cfg_freq = [];
cfg_freq.method = 'wavelet';                % Wavelet-Transformation verwenden
cfg_freq.output = 'powandcsd';              % Power und Kreuzspektraldichte berechnen
cfg_freq.channel = {'F7', 'T7'};            % Die gewünschten Kanäle
cfg_freq.keeptrials = 'yes';                % Trials beibehalten
cfg_freq.foi = 25:1:50;                     % Frequenzen von 25 Hz bis 50 Hz
cfg_freq.toi = -0.4:0.01:0.2;               % Zeitpunkte von -400 ms bis 200 ms
cfg_freq.width = 6;                         % Breite des Morlet-Wavelets
freq_data = ft_freqanalysis(cfg_freq, data); % Frequenzanalyse durchführen

%% WPLI-Analyse
cfg_conn = [];
cfg_conn.method = 'wpli_debiased';          % Debiased wPLI verwenden
cfg_conn.keeptrials = 'yes';                % Trials beibehalten
cfg_conn.channelcmb = {'F7', 'T7'};         % Kanal-Kombination
wpli = ft_connectivityanalysis(cfg_conn, freq_data); % WPLI berechnen

%% Baseline-Zeitfenster definieren (z.B. -400 bis -200 ms)
baseline_time = [-0.4 -0.2];               % Definierte Baseline von -400 ms bis -200 ms
baseline_idx = find(wpli.time >= baseline_time(1) & wpli.time <= baseline_time(2)); % Index für Baseline-Zeiten

%% Z-Wert Berechnung für wPLI-Werte relativ zur Baseline für jede Frequenz
z_wpli = zeros(size(wpli.wpli_debiasedspctrm));  % Initialisiere Z-Wert-Matrix

% Berechnung des Mittelwerts und der Standardabweichung für jede Frequenz separat
for freq_idx = 1:length(wpli.freq)
    % Berechnung des Mittelwerts und der Standardabweichung in der Baseline für die aktuelle Frequenz
    baseline_mean = mean(wpli.wpli_debiasedspctrm(:, freq_idx, baseline_idx), 3);
    baseline_std = std(wpli.wpli_debiasedspctrm(:, freq_idx, baseline_idx), 0, 3);
    
    % Z-Wert Berechnung für jede Zeit für die aktuelle Frequenz
    z_wpli(:, freq_idx, :) = (wpli.wpli_debiasedspctrm(:, freq_idx, :) - baseline_mean) ./ baseline_std;
end

%% Heatmap der Z-transformierten wPLI-Werte
figure;

% Zeit- und Frequenzvektoren
time_vector = wpli.time;  % Zeitpunkte für den Plot
freq_vector = wpli.freq;  % Frequenzen

% Heatmap der Z-Werte der wPLI
imagesc(time_vector, freq_vector, squeeze(mean(z_wpli, 1)));  % Z-transformierte wPLI
axis xy; % Achse umkehren
xlabel('Zeit (s)');
ylabel('Frequenz (Hz)');
title('Z-transformierte wPLI Heatmap über Zeit und Frequenz (Wavelet)');
colorbar; % Farbskala anzeigen
caxis([-3, 3]); % Farbskala für Z-Werte anpassen, optional





