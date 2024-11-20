% pp_subj_sanity_new.m
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
INPATH_ERP = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
INPATH_COH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\subj_sanity'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH_ERP, INPATH_COH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%variables to edit
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj_erp = dir(fullfile(INPATH_ERP, 'P*'));
dircont_subj_coh = dir(fullfile(INPATH_ERP, 'P*'));

%check if both directories contain the same subjects
same_subj_length = length(dircont_subj_erp) == length(dircont_subj_coh); %the script also checks if subjects are the same (not only the number of subjects), see below

%initialize sanity check variables
ok_subj = {};

%setup progress bar
wb = waitbar(0,'starting pp_erp_analysis.m');

for subj = 1:length(dircont_subj_erp) %loop over subjects
    tic;
    %get current ID
    SUBJ_ERP = erase(dircont_subj_erp(subj).name, '_erp_preprocessed.set'); 
    SUBJ_COH = erase(dircont_subj_coh(subj).name, '_coherence_preprocessed.set');
    if ~strcmp(SUBJ_ERP, SUBJ_COH)
        error('subjects do not match')
    else 
        SUBJ = SUBJ_ERP;
    end
    %update progress bar
    waitbar(subj/length(dircont_subj_erp),wb, [SUBJ ' pp_subj_sanity_new.m'])

    %get number of trials for ERP-preprocessed data
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj_erp(subj).name,'filepath',INPATH_ERP);
    %sanity check
    %check if ID matches dataset
    subj_check_erp = strcmp(SUBJ, erase(EEG.setname, '_erp_preprocessed'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen_erp_processed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        EEG.setname = [SUBJ '_' EVENTS{cond} '_erp_preprocessed'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %get number of trials
        switch EVENTS{cond}
            case 'talk'
                erp_num_trials_talk{subj,1} = SUBJ;
                erp_num_trials_talk{subj,2} = size(EEG.data,3);
            case 'listen'
                erp_num_trials_listen{subj,1} = SUBJ;
                erp_num_trials_listen{subj,2} = size(EEG.data,3);
            otherwise
                error('event unknown')
        end
    end

    %get number of trials for coherence-preprocessed data
    %import data
    %start eeglab
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %import dataset
    EEG = pop_loadset('filename',dircont_subj_coh(subj).name,'filepath',INPATH_COH);
    %sanity check
    %check if ID matches dataset
    subj_check_coh = strcmp(SUBJ, erase(EEG.setname, '_coherence_preprocessed'));
    %rename dataset
    EEG.setname = [SUBJ '_talk_listen_coherence_processed'];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    for cond = 1:length(EVENTS) %loop over condition 
        EEG = pop_selectevent( ALLEEG(1), 'latency','-2<=2','type',EVENTS(cond), ...
            'deleteevents','off','deleteepochs','on','invertepochs','off'); %create dataset for each condition
        EEG.setname = [SUBJ '_' EVENTS{cond} '_coherence_preprocessed'];
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        %get number of trials
        switch EVENTS{cond}
            case 'talk'
                coh_num_trials_talk{subj,1} = SUBJ;
                coh_num_trials_talk{subj,2} = size(EEG.data,3);
            case 'listen'
                coh_num_trials_listen{subj,1} = SUBJ;
                coh_num_trials_listen{subj,2} = size(EEG.data,3);
            otherwise
                error('event unknown')
        end
    end
    %sanity checks
    subj_time = toc;
    ok_subj{subj,1} = SUBJ;
    ok_subj{subj,2} = subj_check_erp;
    ok_subj{subj,3} = subj_check_coh;
    ok_subj{subj,4} = subj_time;
end

%create tables
coh_num_trials_talk_table = cell2table(coh_num_trials_talk, 'VariableNames',{'subj','n_coh_talk'});
coh_num_trials_listen_table = cell2table(coh_num_trials_listen, 'VariableNames',{'subj','n_coh_listen'});
erp_num_trials_talk_table = cell2table(erp_num_trials_talk, 'VariableNames',{'subj','n_erp_talk'});
erp_num_trials_listen_table = cell2table(erp_num_trials_listen, 'VariableNames',{'subj','n_erp_listen'});

%create summary table
num_trials_all = join(coh_num_trials_talk_table, coh_num_trials_listen_table);
num_trials_all = join(num_trials_all, erp_num_trials_talk_table);
num_trials_all = join(num_trials_all, erp_num_trials_listen_table);

%sanity check
%see if length of the table matches length of dirconts
check_table = height(num_trials_all) == length(dircont_subj_erp);

%look for subjects with too little trials
ok_num_trials = ~any(num_trials_all{:, 2:end} < 50, 2);
ok_num_trials_subj = num_trials_all.subj(ok_num_trials);

%% Plot number of trials

%extract needed data
subjects = num_trials_all.subj;
conds = {'n_coh_talk', 'n_coh_listen', 'n_erp_talk', 'n_erp_listen'};

% Erstelle die 4 Subplots
figure;

for cond = 1:length(conds)
    values = num_trials_all{:, conds{cond}};
    bar_colors = repmat([0 0 1], length(values), 1); 
    bar_colors(values < 50, :) = repmat([1 0 0], sum(values < 50), 1); 
    subplot(2, 2, cond);
    b = bar(values, 'FaceColor', 'flat'); 
    b.CData = bar_colors; 
    
    xticks(1:length(subjects));
    xticklabels(subjects);
    xlabel('Subjects');
    ylabel('Number of Trials');
    yline(50)
    title(conds{cond}, 'Interpreter', 'none');
end

% Layout optimieren
sgtitle('Number of Trials per Condition');

set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, '_subj_sanity_num_trials.png'))

%% Move clean participants to new folder
%ERP-preprocessed data
subjects = ok_num_trials_subj;

ALLDATAPATH = INPATH_ERP; 
CLEANDATAPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc_clean'); 

%sanity check
%check if folders exist
pp_check_folder_TD(CLEANDATAPATH)
%move current files to archive folder
pp_clean_up_folder_TD(CLEANDATAPATH)

%get all files 
all_files = dir(fullfile(ALLDATAPATH, '*.set')); 

%loop through each file and check for matches
for i = 1:length(all_files)
    file_name = all_files(i).name;
    %check if the file name contains any identifier from the cell array of clean subjects
    for j = 1:length(subjects)
        if contains(file_name, subjects{j})
            %if a match is found, move the file
            copyfile(fullfile(ALLDATAPATH, file_name), fullfile(CLEANDATAPATH, file_name));
        end
    end
end

length_clean_erp = length(dir(fullfile(CLEANDATAPATH, '*.set')));

%coherence-preprocessed data
subjects = ok_num_trials_subj;

ALLDATAPATH = INPATH_COH; 
CLEANDATAPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_coherence_proc_clean'); 

%sanity check
%check if folders exist
pp_check_folder_TD(CLEANDATAPATH)
%move current files to archive folder
pp_clean_up_folder_TD(CLEANDATAPATH)

%get all files 
all_files = dir(fullfile(ALLDATAPATH, '*.set')); 

%loop through each file and check for matches
for i = 1:length(all_files)
    file_name = all_files(i).name;
    %check if the file name contains any identifier from the cell array of clean subjects
    for j = 1:length(subjects)
        if contains(file_name, subjects{j})
            %if a match is found, move the file
            copyfile(fullfile(ALLDATAPATH, file_name), fullfile(CLEANDATAPATH, file_name));
        end
    end
end

length_clean_coh = length(dir(fullfile(CLEANDATAPATH, '*.set')));

%sanity check
%check whether the clean data folders contain as many files as indentified clean subjects
ok_final_dircont_length = length_clean_erp == length_clean_coh && length_clean_erp == length(ok_num_trials_subj);
switch ok_final_dircont_length
    case 1
        disp('all directories OK')
    otherwise
        error('number of files not matching')
end

%% End of processing

%display sanity check variables
ok_subj = cell2table(ok_subj, 'VariableNames',{'subj','subj_check_ID_erp','subj_check_ID_coh','time'})
writetable(ok_subj,fullfile(OUTPATH, '_subj_sanity_ok_subj.xlsx'))

check_done = 'OK'
save(fullfile(OUTPATH, '_subj_sanity_data.mat'), 'check_done', 'check_table', 'ok_num_trials_subj', 'ok_final_dircont_length')

close(wb)
close all


