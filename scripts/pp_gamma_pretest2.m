% pp_gamma_pretest2.m
%
% Pretest with the goal of identifying the presence (or non-presence) of
% gamma (and theta) activity. Using raw data. 
% Only verdical talk condition
% Only for one subject.
%
% DISCONTINUED, see coherence_pretest.m
%
% Tim Dressler, 12.08.2024

clear
close all
clc

%variables to edit
SUBJ = '136'; %122, 125 or 127 available
CHAN = 'Pz';
EVENT = {'S 64'};
EPO_FROM = -0.4;
EPO_TILL = 0.6;
LCF = 5;
HCF = 60;
BL_FROM = -400;
BL_TILL= -200;
THRESH = 75;
SD_PROB = 3;
TF_FROM = -400;
TF_TILL = 500;

%setup paths
MAINPATH = "C:/Users/timdr/OneDrive/Uni_Oldenburg/3_Semester\Module/Pratical_Project/Analysis";
INPATH = fullfile(MAINPATH,['data/raw_data/pp_data_gamma_pretest_raw/P' SUBJ]);
addpath("C:/Users/timdr/OneDrive/Uni_Oldenburg/3_Semester\Module/Pratical_Project/Analysis/neucodis/functions")

%start eeglab & load data
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadbv(INPATH, ['av_P' SUBJ '_C_0001.vhdr'], [], []);

%filter
EEG = pop_eegfiltnew(EEG, 'locutoff',LCF,'hicutoff',HCF);

%epoch
EEG = pop_epoch( EEG,   EVENT  , [EPO_FROM        EPO_TILL], 'epochinfo', 'yes');

%baseline
EEG = pop_rmbase( EEG, [BL_FROM BL_TILL] ,[]);

%treshhold & probability
EEG = pop_eegthresh(EEG,1,[1:64] ,-THRESH,THRESH,-0.2,0.798,2,0);
EEG = pop_jointprob(EEG,1,[1:64] ,SD_PROB,0,0,0,[],0);
EEG = pop_rejkurt(EEG,1,[1:64] ,SD_PROB,0,0,0,[],0);
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);

%find Channel index
chani = find(strcmp({EEG.chanlocs.labels}, CHAN));

%ERSP (Wavelet)
figure; pop_newtimef( EEG, 1, chani, [TF_FROM  TF_TILL], [3         0.8] , 'topovec', chani, ...
    'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', CHAN, 'baseline',[-100 0], ...
    'plotphase', 'off', 'scale', 'abs', 'padratio', 1, 'winsize', 200);

%FastFFT (own function)
%%pp_fastfft_TD(EEG.data, EEG.srate, chani, 'wins', 0.8, 'overlap', 0.99, 'fig', 1, 'timevec', EEG.times, 'ylim_min', 15, 'ylim_max', 50)

%%
%TF-Topoplots

for elec = 1:EEG.nbchan
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = pop_newtimef(EEG, ...
    1, elec, [EEG.xmin EEG.xmax]*1000, [3 0.5], 'maxfreq', 50, 'padratio', 16, ...
    'plotphase', 'off', 'timesout', 60, 'alpha', .05, 'plotersp','off', 'plotitc','off');
    if elec == 1  %create empty arrays if first electrode
        allersp = zeros([ size(ersp) EEG.nbchan]);
        allitc = zeros([ size(itc) EEG.nbchan]);
        allpowbase = zeros([ size(powbase) EEG.nbchan]);
        alltimes = zeros([ size(times) EEG.nbchan]);
        allfreqs = zeros([ size(freqs) EEG.nbchan]);
        allerspboot = zeros([ size(erspboot) EEG.nbchan]);
        allitcboot = zeros([ size(itcboot) EEG.nbchan]);
    end;
    allersp (:,:,elec) = ersp;
    allitc (:,:,elec) = itc;
    allpowbase (:,:,elec) = powbase;
    alltimes (:,:,elec) = times;
    allfreqs (:,:,elec) = freqs;
    allerspboot (:,:,elec) = erspboot;
    allitcboot (:,:,elec) = itcboot;
end;

%Topoplot
figure;
tftopo(allersp,alltimes(:,:,1),allfreqs(:,:,1), ...
    'timefreqs', [50 40; 60 40; 53 42; 100 42], 'chanlocs', EEG.chanlocs, 'showchan', chani)

%%
eeglab rebuild




