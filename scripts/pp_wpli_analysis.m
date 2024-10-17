% pp_wpli_analysis.m
%
% Description
%
% Tim Dressler, 11.09.2024

clear
close all
clc

%selection of electrodes
FRONTAL_R = {'F9', 'F7', 'F5', 'F3'};
TEMPORAL_R = {'T7', 'TP9', 'CP5', 'FC5'};
OCCIPITAL_R = {'O9', 'O1', 'PO9', 'PO7'};

%check if number of electrodes is the same for each lobe
switch length(FRONTAL_R) == length(TEMPORAL_R) && length(FRONTAL_R) == length(OCCIPITAL_R)
    case true
        disp('Electrodes OK')
        NUM_ELE = length(FRONTAL_R);
    otherwise
        error('Electrodes not OK')
end

%initialize electrode pairs variable
PAIRS_L = {};

%generate electrode pairings
r = 1;
for s = 1:NUM_ELE
    for v = 1:NUM_ELE
        PAIRS_L{r, 1} = FRONTAL_R{s};
        PAIRS_L{r, 2} = TEMPORAL_R{v};
        r = r+1;
    end
end

%initialize electrode pairs variable (control)
PAIRS_L_CONTROL = {};

%generate electrode pairings (control)
r = 1;
for s = 1:NUM_ELE
    for v = 1:NUM_ELE
        PAIRS_L_CONTROL{r, 1} = FRONTAL_R{s};
        PAIRS_L_CONTROL{r, 2} = OCCIPITAL_R{v};
        r = r+1;
    end
end

%import data

%convert data to fieldtrip format



%main analysis
cfg_freq = [];
cfg_freq.method = 'mtmfft';
cfg_freq.output = 'powandcsd';
cfg_freq.keeptrials = 'yes'; % do not return an average of all trials for subsequent wpli analysis
cfg_freq.taper = 'dpss';
cfg_freq.tapsmofrq = 1;
cfg_freq.foilim = [0 48];

wpli_ALL = [];

for s = 1:length(PAIRS_L)
    cfg_freq.channel = {PAIRS_L{s,1} PAIRS_L{s,2}};
    freq_data = ft_freqanalysis(cfg_freq, data);
    cfg_conn = [];
    cfg_conn.method = 'wpli_debiased';
    wpli = ft_connectivityanalysis(cfg_conn, freq_data);
    wpli = ft_checkdata(wpli, 'cmbrepresentation', 'full','datatype','freq');
    wpli_ALL(:,:,end+1) = wpli;
end

wpli_MEAN = mean (wpl_ALL,3);






