% pp_subj_sanity.m
%
% Perform subject sanity check.
% Using preprocessed data.
% Preprocessing done with pp_erp_pre_proc.m & pp_coherence_pre_proc.m.
% Check if subjects had enough valid trials for coherence & ERP analysis.
% Stores subjects with enough valid trials in separate folder.
% Subjects were removed from analysis if they had less than 50 trials in
% *either* preprocessing script. Therefore, the subjects analyzed are the
% same for both ERP and coherence analysis.
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
INPATH_COH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
INPATH_ERP = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\subj_sanity'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH_COH, INPATH_ERP, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%variables to edit
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj_COH = dir(fullfile(INPATH_COH, 'P*'));
dircont_subj_ERP= dir(fullfile(INPATH_ERP, 'P*'));

%check if both directories contain the same subjects
same_subj = all([dircont_subj_ERP.name] == [dircont_subj_COH.name]);
same_subj_length = length({dircont_subj_ERP.name}) == length({dircont_subj_COH.name})
switch same_subj && same_subj_length
    case 1
        disp('same subjects present in both directories')
    case 2
        error('stopped processing - different subjects in directories')
end

%initialize sanity check variables
ok_subj = {};

%setup progress bar
wb = waitbar(0,'starting pp_subj_sanity2.m');

r = 1;
for subj = 1:length(dircont_subj_ERP) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj_ERP(subj).name, '_erp_preprocessed.set');
    %update progress bar
    waitbar(subj/length(dircont_subj_ERP),wb, [SUBJ 'ERP pp_subj_sanity.m'])

    %import ERP-preprocessed file
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj_ERP(subj).name,'filepath',INPATH_ERP);
    %sanity check
    %check if ID matches dataset
    subj_check_ERP = strcmp(SUBJ, erase(EEG.setname, '_erp_preprocessed'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        %get number of trials
        num_trials = size(EEG.data,3);
        %store number of trials
        switch EVENTS{cond}
            case 'talk'
                all_num_trials_talk_ERP{subj,1} = SUBJ;
                all_num_trials_talk_ERP{subj,2} = num_trials;
                all_num_trials_talk_ERP{subj,3} = num_trials < 50;
            case 'listen'
                all_num_trials_listen_ERP{subj,1} = SUBJ;
                all_num_trials_listen_ERP{subj,2} = num_trials;
                all_num_trials_listen_ERP{subj,3} = num_trials < 50;
            otherwise
                error('event unknown')
        end
    end

    %import coherence-preprocessed file
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj_COH(subj).name,'filepath',INPATH_COH);
    %sanity check
    %check if ID matches dataset
    subj_check_COH = strcmp(SUBJ, erase(EEG.setname, '_coherence_preprocessed'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        %get number of trials
        num_trials = size(EEG.data,3);
        %store number of trials
        switch EVENTS{cond}
            case 'talk'
                all_num_trials_talk_COH{subj,1} = SUBJ;
                all_num_trials_talk_COH{subj,2} = num_trials;
                all_num_trials_talk_COH{subj,3} = num_trials < 50;
            case 'listen'
                all_num_trials_listen_COH{subj,1} = SUBJ;
                all_num_trials_listen_COH{subj,2} = num_trials;
                all_num_trials_listen_COH{subj,3} = num_trials < 50;
            otherwise
                error('event unknown')
        end
    end

    %sanity checks
    subj_time = toc;
    ok_subj{subj,1} = SUBJ;
    ok_subj{subj,2} = subj_check_ERP;
    ok_subj{subj,3} = subj_check_COH;
    ok_subj{subj,4} = subj_time;
end

%check if any subject has a too little amount of trials
%for ERP-preprocessed data
any_talk_critical_ERP = any([all_num_trials_talk_ERP{:,2}] < 50);
any_listen_critical_ERP = any([all_num_trials_listen_ERP{:,2}] < 50);
%for coherence-preprocessed data
any_talk_critical_COH = any([all_num_trials_talk_COH{:,2}] < 50);
any_listen_critical_COH = any([all_num_trials_listen_COH{:,2}] < 50);

%get table with number of trials
%for ERP-preprocessed data
num_trials_talk_ERP = cell2table(all_num_trials_talk_ERP, 'VariableNames',{'subj','trials_talk', 'critical_talk'});
num_trials_listen_ERP = cell2table(all_num_trials_listen_ERP, 'VariableNames',{'subj','trials_listen', 'critical_listen'});

num_trials_ERP = join(num_trials_talk_ERP, num_trials_listen_ERP);
writetable(num_trials_ERP,fullfile(OUTPATH, '_num_trials_ERP.xlsx'))
%for coherence-preprocessed data
num_trials_talk_COH = cell2table(all_num_trials_talk_COH, 'VariableNames',{'subj','trials_talk', 'critical_talk'});
num_trials_listen_COH = cell2table(all_num_trials_listen_COH, 'VariableNames',{'subj','trials_listen', 'critical_listen'});

num_trials_COH = join(num_trials_talk_COH, num_trials_listen_COH);
writetable(num_trials_COH,fullfile(OUTPATH, '_num_trials_COH.xlsx'))

%plot number of trials
% Plotting number of trials with colored bars
figure;
set(gcf, 'Position', [100, 100, 1200, 800]);

% Extract data for plotting
subjects = num_trials_ERP.subj;
trials_talk_ERP = num_trials_ERP.trials_talk;
trials_listen_ERP = num_trials_ERP.trials_listen;
trials_talk_COH = num_trials_COH.trials_talk;
trials_listen_COH = num_trials_COH.trials_listen;

% Define threshold
threshold = 50;

% Subplot 1: ERP - Talk
subplot(2, 2, 1);
bar(trials_talk_ERP, 'FaceColor', 'flat');
for i = 1:length(trials_talk_ERP)
    if trials_talk_ERP(i) < threshold
        color = [1, 0, 0]; % Red for < 50
    else
        color = [0, 0.4470, 0.7410]; % Default MATLAB blue
    end
    set(get(gca, 'Children'), 'CData', i, 'FaceColor', color);
end
title('ERP - Talk');
xlabel('Subjects');
ylabel('Number of Trials');
xticks(1:length(subjects));
xticklabels(subjects);
xtickangle(45);

% Subplot 2: ERP - Listen
subplot(2, 2, 2);
bar(trials_listen_ERP, 'FaceColor', 'flat');
for i = 1:length(trials_listen_ERP)
    if trials_listen_ERP(i) < threshold
        color = [1, 0, 0]; % Red for < 50
    else
        color = [0, 0.4470, 0.7410]; % Default MATLAB blue
    end
    set(get(gca, 'Children'), 'CData', i, 'FaceColor', color);
end
title('ERP - Listen');
xlabel('Subjects');
ylabel('Number of Trials');
xticks(1:length(subjects));
xticklabels(subjects);
xtickangle(45);

% Subplot 3: Coherence - Talk
subplot(2, 2, 3);
bar(trials_talk_COH, 'FaceColor', 'flat');
for i = 1:length(trials_talk_COH)
    if trials_talk_COH(i) < threshold
        color = [1, 0, 0]; % Red for < 50
    else
        color = [0, 0.4470, 0.7410]; % Default MATLAB blue
    end
    set(get(gca, 'Children'), 'CData', i, 'FaceColor', color);
end
title('Coherence - Talk');
xlabel('Subjects');
ylabel('Number of Trials');
xticks(1:length(subjects));
xticklabels(subjects);
xtickangle(45);

% Subplot 4: Coherence - Listen
subplot(2, 2, 4);
bar(trials_listen_COH, 'FaceColor', 'flat');
for i = 1:length(trials_listen_COH)
    if trials_listen_COH(i) < threshold
        color = [1, 0, 0]; % Red for < 50
    else
        color = [0, 0.4470, 0.7410]; % Default MATLAB blue
    end
    set(get(gca, 'Children'), 'CData', i, 'FaceColor', color);
end
title('Coherence - Listen');
xlabel('Subjects');
ylabel('Number of Trials');
xticks(1:length(subjects));
xticklabels(subjects);
xtickangle(45);

% Overall layout
sgtitle('Number of Trials by Subject and Condition');



%% End of processing

%display sanity check variables
ok_subj = cell2table(ok_subj, 'VariableNames',{'subj','subj_check_ID_ERP','subj_check_ID_COH','time'})
writetable(ok_subj,fullfile(OUTPATH, '_subj_sanity_ok_subj.xlsx'));

check_done = 'OK'
%%save(fullfile(OUTPATH, '_subj_sanity_data.mat'), 'check_done', 'any_listen_critical', 'any_talk_critical')

close(wb)