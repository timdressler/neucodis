% pp_subj_sanity.m
%
% Perform subject sanity check.
% Using preprocessed data.
% Preprocessing done with pp_erp_pre_proc.m & pp_coherence_pre_proc.m.
% Check if subjects had enough valid trials for coherence & ERP analysis.
%
% Tim Dressler, 18.11.2024

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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\subj_sanity'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%variables to edit
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
ok_subj = {};

%setup progress bar
wb = waitbar(0,'starting pp_subj_sanity.m');

r = 1;
for subj = 1:length(dircont_subj) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_erp_preprocessed.set');
    %update progress bar
    waitbar(subj/length(dircont_subj),wb, [SUBJ ' pp_subj_sanity.m'])
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj(subj).name,'filepath',INPATH);
    %sanity check
    %check if ID matches dataset
    subj_check = strcmp(SUBJ, erase(EEG.setname, '_erp_preprocessed'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        EEG.setname = [SUBJ EVENTS{cond}];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %get number of trials
        num_trials = size(EEG.data,3);
        %store number of trials
        switch EVENTS{cond}
            case 'talk'
                all_num_trials_talk{subj,1} = SUBJ;
                all_num_trials_talk{subj,2} = num_trials;
                all_num_trials_talk{subj,3} = num_trials < 50;
            case 'listen'
                all_num_trials_listen{subj,1} = SUBJ;
                all_num_trials_listen{subj,2} = num_trials;
                all_num_trials_listen{subj,3} = num_trials < 50;
            otherwise
                error('event unknown')
        end
    end
    %sanity checks
    subj_time = toc;
    ok_subj{subj,1} = SUBJ;
    ok_subj{subj,2} = subj_check;
    ok_subj{subj,3} = subj_time;
end

%check if any subject has a too little amount of trials
any_talk_critical = any([all_num_trials_talk{:,2}] < 50);
any_listen_critical = any([all_num_trials_listen{:,2}] < 50);

num_trials_talk = cell2table(all_num_trials_talk, 'VariableNames',{'subj','trials_talk', 'critical_talk'});
num_trials_listen = cell2table(all_num_trials_listen, 'VariableNames',{'subj','trials_listen', 'critical_listen'});

num_trials = join(num_trials_talk, num_trials_listen);
writetable(num_trials,fullfile(OUTPATH, '_num_trials.xlsx'))

%plot number of trials
%talk condition
figure;
subplot(121)
bar([all_num_trials_talk{:,2}])
hold on
yline(50)
hold off
title('talk')

%listen condition
subplot(122)
bar([all_num_trials_listen{:,2}])
hold on
yline(50)
hold off
title('listen')

sgtitle('Numer of trials for all subjects')

set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'num_trials_listen.png'))

%% End of processing

%display sanity check variables
ok_subj = cell2table(ok_subj, 'VariableNames',{'subj','subj_check_ID','time'})
writetable(ok_subj,fullfile(OUTPATH, '_subj_sanity_ok_subj.xlsx'));

check_done = 'OK'
save(fullfile(OUTPATH, '_subj_sanity_data.mat'), 'check_done', 'any_listen_critical', 'any_talk_critical')

close(wb)

