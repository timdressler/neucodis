% pp_erp_analysis.m
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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
CHAN = 'Cz';
ERP_FROM = 50;
ERP_TILL = 150;
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
marked_subj = {};
ok_subj = {};

r = 1;
for subj = 1:length(dircont_subj) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_erp_preprocessed.set');
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
    %sanity check
    %check if ID matches dataset
    subj_check = strcmp(SUBJ, erase(EEG.setname, '_erp_preprocessed'));
    %get channel ID
    chani = find(strcmp({EEG.chanlocs.labels}, CHAN));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        EEG.setname = [SUBJ EVENTS{cond}];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %get ERP
        ERP = mean(EEG.data, 3);
        %setup ERP analysis
        [~,win_start] = min(abs(EEG.times-ERP_FROM));
        [~,win_end] = min(abs(EEG.times-ERP_TILL));
        [~,t_zero] = min(abs(EEG.times));
        %get N100 amplitude
        maxERPamp = min(ERP(chani,win_start:win_end));
        %get N100 latency
        maxERPsam = find(ERP(chani,:) == maxERPamp);
        maxERPlat = EEG.times(maxERPsam);
        %store ERPs
        switch EVENTS{cond}
            case 'talk'
                all_ERP_talk(:,:, subj) = ERP;
            case 'listen'
                all_ERP_listen(:,:, subj) = ERP;
            otherwise
                error('event unknown')
        end
        %store variables
        all_ERP{r,1} = ERP_FROM;
        all_ERP{r,2} = ERP_TILL;
        all_ERP{r,3} = maxERPamp;
        all_ERP{r,4} = maxERPlat;
        all_ERP{r,5} = EVENTS{cond};
        all_ERP{r,6} = SUBJ;
        r = r+1;
    end
    %sanity checks
    subj_time = toc;
    ok_subj{subj,1} = SUBJ;
    ok_subj{subj,2} = subj_check;
    ok_subj{subj,3} = subj_time;
end

%get grandaverage ERP


%store variable in table
all_ERP_table = table(all_ERP(:,6),all_ERP(:,1), all_ERP(:,2), all_ERP(:,3), ...
    all_ERP(:,4), all_ERP(:,5), 'VariableNames',{'subj','erp_from', 'erp_till', 'erp_amp', 'erp_lat', 'cond'})


ok_subj


plot(EEG.times, all_ERP_listen(18,:,1))
hold on
plot(EEG.times, all_ERP_talk(18,:,1))
legend
hold off 



















