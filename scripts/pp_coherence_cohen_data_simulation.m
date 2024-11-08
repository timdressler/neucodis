% pp_wPLI_analysis.m
%
% Data simulation for sanity-check of coherence analysis
% Data can be only noise or contain simulated coherence
%
% Originally by Dr. Micheal X Cohen
% Cohen, MX (2014). Effects of time lag and frequency matching on
% phase-based connectivity. Journal of Neuroscience Methods.
%
% Adapted from Tim Dressler, 19.09.2024

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
OUTPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_simulated\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)

%variables to edit
SIM_COH = 1; %1 for simulated coherence (only in talk condition) %2 for simulated coherence (in both conditions) %0 for noise data
COH_START = -0.2; %in s %only has an effect if SIM_COH = 1

%load participant and extract number of trials to create simulated dataset equivalent to the original data
%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
ok_subj = {};

switch SIM_COH
    case 0 %noise data
        for subj = 1:length(dircont_subj) %loop over subjects
            %general setup
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
            %extract neeeded information
            %time vector
            time = EEG.times;
            %sample rate
            srate = EEG.srate;
            %number of trials
            ntrials = size(EEG.data,3);
            %event data
            event = EEG.event;
            %epoch data
            epoch = EEG.epoch;

            %simulation setup
            %center frequency in Hz
            centfreq = 35;
            %mean phase difference between two dipoles (state in ms, converted later to radians)
            phasedif_ms = 5; % 5 or 25
            %simulationType: 1 (equal stationary frequencies)
            %                2 (unequal stationary frequencies)
            %                3 (nonstationary frequencies)
            simulationType = 1;
            %dispersion of frequency nonstationarity, in Hz. Values between 1 and 5 are reasonable.
            %it has an effect only when stimulationType is set to 3.
            freqDist = 2;
            %number of trials
            ntrials = ntrials;
            %initial setup and load in necessary files
            %mat file containing leadfield and channel locations
            load lfchans
            %sampling rate
            srate = srate;
            %time for simulation (in seconds)
            time  = time;
            phasedif = 2*pi*centfreq*phasedif_ms/1000;
            ntime  = length(time);
            nchans = length(chanlocs);
            %indices of dipole locations (probably best not to change these)
            dipoleOCC =   94;
            dipolePFC = 1720;
            %use X, Y, or Z oriented dipole (1, 2, or 3, respectively).
            %in the paper, Z was used.
            whichOrientation = 3;
            %to see the scalp dipole projections, use the following code
            %topoplot(squeeze(lf.Gain(:,whichOrientation,dipolePFC)),chanlocs);
            %intialize output matrices
            simulatedEEG = zeros(nchans,ntime,ntrials); %simulated electrode data
            sourceTimeSeries = zeros(ntime,2,ntrials); %data at source dipoles

            %start simulation
            for triali=1:ntrials
                %data comprise all noise
                data = randn(ntime,size(lf.Gain,3))./20;
                %simulated EEG data
                simulatedEEG(:,:,triali) = (data*squeeze(lf.Gain(:,whichOrientation,:))')';
                %get actual source time series
                sourceTimeSeries(:,:,triali) = data(:,[dipolePFC dipoleOCC]);
            end
            %also compute laplacian
            % % simulatedLap = laplacian_perrinX(simulatedEEG,[chanlocs.X],[chanlocs.Y],[chanlocs.Z]);
            %average reref (leadfield assumes average reference)
            simulatedEEG = bsxfun(@minus,simulatedEEG,mean(simulatedEEG,1));
            %convert simulated data to eeglab format
            EEG.data = simulatedEEG;
            EEG.times = time;
            EEG.srate = srate;
            EEG.chanlocs = chanlocs;
            EEG.event = event;
            EEG.epoch = epoch;
            EEG.nbchan = size(EEG.data, 1);
            EEG.trials = ntrials;
            EEG.pnts = size(EEG.data, 2);
            EEG.sim = 'Noise';
            %store dataset
            EEG.setname = [SUBJ '_simulated_talk'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            %end of simulation

            %save dataset
            EEG.setname = [SUBJ '_simulated'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            EEG = pop_saveset(EEG, 'filename',[SUBJ '_simulated.set'],'filepath', OUTPATH);

            %sanity checks
            subj_time = toc;
            ok_subj{subj,1} = [SUBJ '_simulated'];
            ok_subj{subj,2} = subj_check;
            ok_subj{subj,3} = subj_time;
        end

    case 1 %simulated coherence (only in talk condition)
        for subj = 1:length(dircont_subj) %loop over subjects
            %general setup
            tic;
            %get current ID
            SUBJ = erase(dircont_subj(subj).name, '_coherence_preprocessed.set');

            %simulate talk trials (with coherence)
            %import data
            %start eeglab
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            %import dataset
            EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
            %sanity check
            %check if ID matches dataset
            subj_check = strcmp(SUBJ, erase(EEG.setname, '_coherence_preprocessed'));
            %extract neeeded information
            %time vector
            time = EEG.times;
            %sample rate
            srate = EEG.srate;
            %number of trials
            ntrials = size(EEG.data,3);
            %event data
            event = EEG.event;
            %epoch data
            epoch = EEG.epoch;

            %simulation setup
            %center frequency in Hz
            centfreq = 35;
            %mean phase difference between two dipoles (state in ms, converted later to radians)
            phasedif_ms = 5; % 5 or 25
            %simulationType: 1 (equal stationary frequencies)
            %                2 (unequal stationary frequencies)
            %                3 (nonstationary frequencies)
            simulationType = 1;
            %dispersion of frequency nonstationarity, in Hz. Values between 1 and 5 are reasonable.
            %it has an effect only when stimulationType is set to 3.
            freqDist = 2;
            %number of trials
            ntrials = ntrials;
            %initial setup and load in necessary files
            %mat file containing leadfield and channel locations
            load lfchans
            %sampling rate
            srate = srate;
            %time for simulation (in seconds)
            time  = time;
            phasedif = 2*pi*centfreq*phasedif_ms/1000;
            ntime  = length(time);
            nchans = length(chanlocs);
            %indices of dipole locations (probably best not to change these)
            dipoleOCC =   94;
            dipolePFC = 1720;
            %use X, Y, or Z oriented dipole (1, 2, or 3, respectively).
            %in the paper, Z was used.
            whichOrientation = 3;
            %to see the scalp dipole projections, use the following code
            %topoplot(squeeze(lf.Gain(:,whichOrientation,dipolePFC)),chanlocs);
            %intialize output matrices
            simulatedEEG = zeros(nchans,ntime,ntrials); %simulated electrode data
            sourceTimeSeries = zeros(ntime,2,ntrials); %data at source dipoles

            %start simulation
            for triali=1:ntrials
                %data comprise all noise
                data = randn(ntime,size(lf.Gain,3))./20;

                switch simulationType
                    case 1
                        %case 1: both dipoles oscillate at centfreq Hz
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq;
                        [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
                    case 2
                        %case 2: one dipole oscillates faster; both oscillators are stationary
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq*1.2;
                        [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
                    case 3
                        %case 3: Both oscillators are nonstationary
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq;
                        %time-varying frequency modulation
                        freqmod1 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod1(ntime+1:end) = [];
                        freqmod2 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod2(ntime+1:end) = [];
                        %frequency coefficients for generating arbitrary 'chirp' signal
                        k1 = (centfreq/srate)*2*pi/trialFreq1;
                        k2 = (centfreq/srate)*2*pi/trialFreq2;
                end
                %gaussian window used for tapering sine waves
                gausWindow = exp( (-((time/1000)-COH_START).^2)/.01 ); %WATCHOUT %.01 was originally .1, thought to control the time range of coherence
                %create signals
                data(:,dipolePFC) = data(:,dipolePFC)' + sin(2*pi.*trialFreq1.*(time/1000) + k1*cumsum(freqmod1) + rand*.1*pi         ) .* gausWindow;
                tempts            = data(:,dipolePFC)' + sin(2*pi.*trialFreq2.*(time/1000) + k2*cumsum(freqmod2) + rand*.1*pi+phasedif) .* gausWindow;
                data(:,dipoleOCC) = data(:,dipoleOCC)' + tempts.*gausWindow;

                %simulated EEG data
                simulatedEEG(:,:,triali) = (data*squeeze(lf.Gain(:,whichOrientation,:))')';
                %get actual source time series
                sourceTimeSeries(:,:,triali) = data(:,[dipolePFC dipoleOCC]);
            end
            %also compute laplacian
            % % simulatedLap = laplacian_perrinX(simulatedEEG,[chanlocs.X],[chanlocs.Y],[chanlocs.Z]);
            %average reref (leadfield assumes average reference)
            simulatedEEG = bsxfun(@minus,simulatedEEG,mean(simulatedEEG,1));
            %convert simulated data to eeglab format
            EEG.data = simulatedEEG;
            EEG.times = time;
            EEG.srate = srate;
            EEG.chanlocs = chanlocs;
            EEG.event = event;
            EEG.epoch = epoch;
            EEG.nbchan = size(EEG.data, 1);
            EEG.trials = ntrials;
            EEG.pnts = size(EEG.data, 2);
            EEG.sim = 'Coh_Talk';
            %select only listen trials
            EEG = pop_selectevent( EEG, 'latency','-2<=2','type',{'talk'},'deleteevents','off','deleteepochs','on','invertepochs','off');
            %store dataset
            EEG.setname = [SUBJ '_simulated_talk'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

            %simulate listen trials (without coherence)
            %import dataset
            EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
            %extract neeeded information
            %time vector
            time = EEG.times;
            %sample rate
            srate = EEG.srate;
            %number of trials
            ntrials = size(EEG.data,3);
            %event data
            event = EEG.event;
            %epoch data
            epoch = EEG.epoch;

            %simulation setup
            %center frequency in Hz
            centfreq = 35;
            %mean phase difference between two dipoles (state in ms, converted later to radians)
            phasedif_ms = 5; % 5 or 25
            %simulationType: 1 (equal stationary frequencies)
            %                2 (unequal stationary frequencies)
            %                3 (nonstationary frequencies)
            simulationType = 1;
            %dispersion of frequency nonstationarity, in Hz. Values between 1 and 5 are reasonable.
            %it has an effect only when stimulationType is set to 3.
            freqDist = 2;
            %number of trials
            ntrials = ntrials;
            %initial setup and load in necessary files
            %mat file containing leadfield and channel locations
            load lfchans
            %sampling rate
            srate = srate;
            %time for simulation (in seconds)
            time  = time;
            phasedif = 2*pi*centfreq*phasedif_ms/1000;
            ntime  = length(time);
            nchans = length(chanlocs);
            %indices of dipole locations (probably best not to change these)
            dipoleOCC =   94;
            dipolePFC = 1720;
            %use X, Y, or Z oriented dipole (1, 2, or 3, respectively).
            %in the paper, Z was used.
            whichOrientation = 3;
            %to see the scalp dipole projections, use the following code
            %topoplot(squeeze(lf.Gain(:,whichOrientation,dipolePFC)),chanlocs);
            %intialize output matrices
            simulatedEEG = zeros(nchans,ntime,ntrials); %simulated electrode data
            sourceTimeSeries = zeros(ntime,2,ntrials); %data at source dipoles

            %start simulation
            for triali=1:ntrials
                %data comprise all noise
                data = randn(ntime,size(lf.Gain,3))./20;
                %simulated EEG data
                simulatedEEG(:,:,triali) = (data*squeeze(lf.Gain(:,whichOrientation,:))')';
                %get actual source time series
                sourceTimeSeries(:,:,triali) = data(:,[dipolePFC dipoleOCC]);
            end
            %also compute laplacian
            % % simulatedLap = laplacian_perrinX(simulatedEEG,[chanlocs.X],[chanlocs.Y],[chanlocs.Z]);
            %average reref (leadfield assumes average reference)
            simulatedEEG = bsxfun(@minus,simulatedEEG,mean(simulatedEEG,1));
            %convert simulated data to eeglab format
            EEG.data = simulatedEEG;
            EEG.times = time;
            EEG.srate = srate;
            EEG.chanlocs = chanlocs;
            EEG.event = event;
            EEG.epoch = epoch;
            EEG.nbchan = size(EEG.data, 1);
            EEG.trials = ntrials;
            EEG.pnts = size(EEG.data, 2);
            EEG.sim = 'Coh_talk';
            %select only listen trials
            EEG = pop_selectevent( EEG, 'latency','-2<=2','type',{'listen'},'deleteevents','off','deleteepochs','on','invertepochs','off');
            %store dataset
            EEG.setname = [SUBJ '_simulated_listen'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            %end of simulation

            %merge dataset
            EEG = pop_mergeset( ALLEEG, [1  2], 0);

            %save dataset
            EEG.setname = [SUBJ '_simulated'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            EEG = pop_saveset(EEG, 'filename',[SUBJ '_simulated.set'],'filepath', OUTPATH);

            %sanity checks
            subj_time = toc;
            ok_subj{subj,1} = [SUBJ '_simulated'];
            ok_subj{subj,2} = subj_check;
            ok_subj{subj,3} = subj_time;
        end
    case 2 %simulated coherence (in both conditions)

        for subj = 1:length(dircont_subj) %loop over subjects
            %general setup
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
            %extract neeeded information
            %time vector
            time = EEG.times;
            %sample rate
            srate = EEG.srate;
            %number of trials
            ntrials = size(EEG.data,3);
            %event data
            event = EEG.event;
            %epoch data
            epoch = EEG.epoch;

            %simulation setup
            %center frequency in Hz
            centfreq = 35;
            %mean phase difference between two dipoles (state in ms, converted later to radians)
            phasedif_ms = 5; % 5 or 25
            %simulationType: 1 (equal stationary frequencies)
            %                2 (unequal stationary frequencies)
            %                3 (nonstationary frequencies)
            simulationType = 1;
            %dispersion of frequency nonstationarity, in Hz. Values between 1 and 5 are reasonable.
            %it has an effect only when stimulationType is set to 3.
            freqDist = 2;
            %number of trials
            ntrials = ntrials;
            %initial setup and load in necessary files
            %mat file containing leadfield and channel locations
            load lfchans
            %sampling rate
            srate = srate;
            %time for simulation (in seconds)
            time  = time;
            phasedif = 2*pi*centfreq*phasedif_ms/1000;
            ntime  = length(time);
            nchans = length(chanlocs);
            %indices of dipole locations (probably best not to change these)
            dipoleOCC =   94;
            dipolePFC = 1720;
            %use X, Y, or Z oriented dipole (1, 2, or 3, respectively).
            %in the paper, Z was used.
            whichOrientation = 3;
            %to see the scalp dipole projections, use the following code
            %topoplot(squeeze(lf.Gain(:,whichOrientation,dipolePFC)),chanlocs);
            %intialize output matrices
            simulatedEEG = zeros(nchans,ntime,ntrials); %simulated electrode data
            sourceTimeSeries = zeros(ntime,2,ntrials); %data at source dipoles

            %start simulation
            for triali=1:ntrials
                %data comprise all noise
                data = randn(ntime,size(lf.Gain,3))./20;
                switch simulationType
                    case 1
                        %case 1: both dipoles oscillate at centfreq Hz
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq;
                        [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
                    case 2
                        %case 2: one dipole oscillates faster; both oscillators are stationary
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq*1.2;
                        [freqmod1,freqmod2,k1,k2] = deal(0); % see case 3
                    case 3
                        %case 3: Both oscillators are nonstationary
                        trialFreq1 = centfreq;
                        trialFreq2 = centfreq;
                        %time-varying frequency modulation
                        freqmod1 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod1(ntime+1:end) = [];
                        freqmod2 = interp(ceil(freqDist*rand(1,round(ntime/40)))-.5-(freqDist/2),41); freqmod2(ntime+1:end) = [];
                        %frequency coefficients for generating arbitrary 'chirp' signal
                        k1 = (centfreq/srate)*2*pi/trialFreq1;
                        k2 = (centfreq/srate)*2*pi/trialFreq2;
                end
                %gaussian window used for tapering sine waves
                gausWindow = exp( (-((time/1000)-COH_START).^2)/.01 ); %WATCHOUT %.01 was originally .1, thought to control the time range of coherence
                %create signals
                data(:,dipolePFC) = data(:,dipolePFC)' + sin(2*pi.*trialFreq1.*(time/1000) + k1*cumsum(freqmod1) + rand*.1*pi         ) .* gausWindow;
                tempts            = data(:,dipolePFC)' + sin(2*pi.*trialFreq2.*(time/1000) + k2*cumsum(freqmod2) + rand*.1*pi+phasedif) .* gausWindow;
                data(:,dipoleOCC) = data(:,dipoleOCC)' + tempts.*gausWindow;
                %simulated EEG data
                simulatedEEG(:,:,triali) = (data*squeeze(lf.Gain(:,whichOrientation,:))')';
                %get actual source time series
                sourceTimeSeries(:,:,triali) = data(:,[dipolePFC dipoleOCC]);
            end
            %also compute laplacian
            % % simulatedLap = laplacian_perrinX(simulatedEEG,[chanlocs.X],[chanlocs.Y],[chanlocs.Z]);
            %average reref (leadfield assumes average reference)
            simulatedEEG = bsxfun(@minus,simulatedEEG,mean(simulatedEEG,1));
            %convert simulated data to eeglab format
            EEG.data = simulatedEEG;
            EEG.times = time;
            EEG.srate = srate;
            EEG.chanlocs = chanlocs;
            EEG.event = event;
            EEG.epoch = epoch;
            EEG.nbchan = size(EEG.data, 1);
            EEG.trials = ntrials;
            EEG.pnts = size(EEG.data, 2);
            EEG.sim = 'Coh_Both';
            %end of simulation

            %save dataset
            EEG.setname = [SUBJ '_simulated'];
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
            EEG = pop_saveset(EEG, 'filename',[SUBJ '_simulated.set'],'filepath', OUTPATH);

            %sanity checks
            subj_time = toc;
            ok_subj{subj,1} = [SUBJ '_simulated'];
            ok_subj{subj,2} = subj_check;
            ok_subj{subj,3} = subj_time;
        end
    otherwise
        error('error')
end

close all

%display sanity check variables
ok_subj
check_done = 'OK'
