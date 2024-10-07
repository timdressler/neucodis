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

%preprocessing dataset C_0001
%load raw dataset (C_0001)
EEG = pop_loadbv(fullfile(INPATH, 'P136'), 'av_P136_C_0001.vhdr', [], []);
%remove not needed channels
EEG = pop_select( EEG, 'rmchannel',{'heogl','heogr','ml','mr', 'Lip'});
%rename events
for s = 2:length(EEG.event)
    if s < length(EEG.event)
        if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64') && ~strcmp(EEG.event(s+1).type, 'S 64')
            EEG.event(s).type = 'talk';
        else
            continue
        end
    elseif s == length(EEG.event)
        if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64')
            EEG.event(s).type = 'talk';
        else
            continue
        end
    end
end
%1-60Hz bandpass
EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',60,'plotfreqz',0);
%run PREP pipeline
%detrend
[EEG, ~] = removeTrend(EEG);
%line noise removal
lineNoiseIn = getPrepDefaults(EEG, 'linenoise');

lineNoiseIn.fPassBand = lineNoiseIn.fPassBand.value;
lineNoiseIn.Fs = lineNoiseIn.Fs.value;
lineNoiseIn.fScanBandWidth = lineNoiseIn.fScanBandWidth.value;
lineNoiseIn.lineFrequencies = 50:50:EEG.srate/2;
lineNoiseIn.lineNoiseChannels = lineNoiseIn.lineNoiseChannels.value;
lineNoiseIn.lineNoiseMethod = lineNoiseIn.lineNoiseMethod.value;
lineNoiseIn.maximumIterations = lineNoiseIn.maximumIterations.value;
lineNoiseIn.p = lineNoiseIn.p.value;
lineNoiseIn.pad = lineNoiseIn.pad.value;
lineNoiseIn.taperBandWidth = lineNoiseIn.taperBandWidth.value;
lineNoiseIn.taperWindowSize = lineNoiseIn.taperWindowSize.value;
lineNoiseIn.taperWindowStep = lineNoiseIn.taperWindowStep.value;
lineNoiseIn.tau = lineNoiseIn.tau.value;

[EEG, lineNoiseOut] = cleanLineNoise(EEG, lineNoiseIn);
%robust rereferencing
[EEG, referenceOut] = performReference(EEG);
%store dataset C_0001
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

%preprocessing dataset C_0005
%load raw dataset (C_0005)
EEG = pop_loadbv(fullfile(INPATH, 'P136'), 'av_P136_C_0001.vhdr', [], []);
%remove not needed channels
EEG = pop_select( EEG, 'rmchannel',{'heogl','heogr','ml','mr', 'Lip'});
%rename events
for s = 2:length(EEG.event)
    if s < length(EEG.event)
        if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64') && ~strcmp(EEG.event(s+1).type, 'S 64')
            EEG.event(s).type = 'listen';
        else
            continue
        end
    elseif s == length(EEG.event)
        if strcmp(EEG.event(s).type, 'S 64') && ~strcmp(EEG.event(s-1).type, 'S 64')
            EEG.event(s).type = 'listen';
        else
            continue
        end
    end
end
%1-60Hz bandpass
EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',60,'plotfreqz',0);
%run PREP pipeline
%detrend
[EEG, ~] = removeTrend(EEG);
%line noise removal
lineNoiseIn = getPrepDefaults(EEG, 'linenoise');

lineNoiseIn.fPassBand = lineNoiseIn.fPassBand.value;
lineNoiseIn.Fs = lineNoiseIn.Fs.value;
lineNoiseIn.fScanBandWidth = lineNoiseIn.fScanBandWidth.value;
lineNoiseIn.lineFrequencies = 50:50:EEG.srate/2;
lineNoiseIn.lineNoiseChannels = lineNoiseIn.lineNoiseChannels.value;
lineNoiseIn.lineNoiseMethod = lineNoiseIn.lineNoiseMethod.value;
lineNoiseIn.maximumIterations = lineNoiseIn.maximumIterations.value;
lineNoiseIn.p = lineNoiseIn.p.value;
lineNoiseIn.pad = lineNoiseIn.pad.value;
lineNoiseIn.taperBandWidth = lineNoiseIn.taperBandWidth.value;
lineNoiseIn.taperWindowSize = lineNoiseIn.taperWindowSize.value;
lineNoiseIn.taperWindowStep = lineNoiseIn.taperWindowStep.value;
lineNoiseIn.tau = lineNoiseIn.tau.value;

[EEG, lineNoiseOut] = cleanLineNoise(EEG, lineNoiseIn);
%robust rereferencing
[EEG, referenceOut] = performReference(EEG);
%store dataset C_0005
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

%preprocessing on merged datasets
%merge datasets
EEG = pop_mergeset( ALLEEG, [1 2], 0);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%run ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
%label ICA components with IC Label Plugin
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, [0 0.2;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
EEG = pop_subcomp( EEG, [], 0);
%compute surface laplacian
EEG.data=laplacian_perrinX(EEG.data,EEG.chanlocs.X,EEG.chanlocs.Y,EEG.chanlocs.Z);
%20-55Hz bandpass
EEG = pop_eegfiltnew(EEG, 'locutoff',20,'hicutoff',55,'plotfreqz',0);
%epoching
EEG = pop_epoch( EEG, {  'talk'  }, [-0.8         0.7], 'epochinfo', 'yes');
%baseline removal
EEG = pop_rmbase( EEG, [-800 -400] ,[]);
%threshold removal
EEG = pop_eegthresh(EEG,1,[1:60] ,-50,50,-0.8,0.699,2,1); %CHECK %changed from +/- 100 to +/- 50
%probability-based removal
EEG = pop_jointprob(EEG,1,[1:60] ,3,3,0,1,0,[],0);
EEG = pop_rejkurt(EEG,1,[1:60] ,3,3,0,1,0,[],0);

%end of preprocessiing

eeglab redraw





