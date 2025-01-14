% pp_erp_analysis.m
%
% Perform ERP analysis.
% Using preprocessed data.
% Preprocessing done with pp_erp_pre_proc.m.
% Using talk & listen conditions.
% Plots grand average ERP over all subjects.
% Exports ERP properties for further analysis in R (see pp_erp_analysis_R.R)
%
% Tim Dressler, 11.09.2024

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
INPATH = fullfile(MAINPATH, 'data\proc_data\pp_data_erp_proc_clean\'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\analysis_data\erp_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%variables to edit
CHAN = 'Cz'; %set to 'Cz' in original analysis
ERP_FROM = 50; %set to 50 in original analysis
ERP_TILL = 150; %set to 150 in original analysis
EVENTS = {'talk', 'listen'};

%get directory content
dircont_subj = dir(fullfile(INPATH, 'P*'));

%initialize sanity check variables
ok_subj = {};

%setup progress bar
wb = waitbar(0,'starting pp_erp_analysis.m');

r = 1;
for subj = 1:length(dircont_subj) %loop over subjects
    tic;
    %get current ID
    SUBJ = erase(dircont_subj(subj).name, '_erp_preprocessed.set');
    %update progress bar
    waitbar(subj/length(dircont_subj),wb, [SUBJ ' pp_erp_analysis.m'])
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
grandaverage_ERP_talk = mean(all_ERP_talk,3);
grandaverage_ERP_listen = mean(all_ERP_listen,3);


%store variable in table
all_ERP_table = table(all_ERP(:,6),all_ERP(:,1), all_ERP(:,2), all_ERP(:,3), ...
    all_ERP(:,4), all_ERP(:,5), 'VariableNames',{'subj','erp_from', 'erp_till', 'erp_amp', 'erp_lat', 'cond'})

%export table
writetable(all_ERP_table,fullfile(OUTPATH, 'erp_analysis.csv'))

%% Plot ERP

%set colors 
main_blue = '#004F9F';
main_blue = sscanf(main_blue(2:end),'%2x%2x%2x',[1 3])/255;
main_red = '#D53D0E';
main_red = sscanf(main_red(2:end),'%2x%2x%2x',[1 3])/255;
main_green = '#00786B';
main_green = sscanf(main_green(2:end),'%2x%2x%2x',[1 3])/255;
light_blue = '#5BC5F2';
light_blue = sscanf(light_blue(2:end),'%2x%2x%2x',[1 3])/255;
main_yellow = '#FDC300';
main_yellow = sscanf(main_yellow(2:end),'%2x%2x%2x',[1 3])/255;

%get values for dynamic plot limits
y_lim_upper = min([min(all_ERP_talk(chani,:,:), [], 'all') min(all_ERP_listen(chani,:,:), [], 'all')])-0.4;
y_lim_lower = max([max(all_ERP_talk(chani,:,:), [], 'all') max(all_ERP_listen(chani,:,:), [], 'all')])+0.4;


%plot
close all
figure;
p(1:1+size(all_ERP_listen,3)-1) = ...
    plot(EEG.times, squeeze(all_ERP_talk(18,:,:)).', ...
    'linestyle','--','color',main_red,'linewidth',0.01);
hold on
p(1+size(all_ERP_listen,3): 1+size(all_ERP_listen,3)+size(all_ERP_talk,3)-1) = ...
    plot(EEG.times, squeeze(all_ERP_listen(18,:,:)).', ...
    'linestyle','--','color',main_blue,'linewidth',0.01);
p(1+size(all_ERP_listen,3)+size(all_ERP_talk,3)) = plot(EEG.times, grandaverage_ERP_talk(chani,:),'color', main_red, 'LineWidth',2);
p(1+size(all_ERP_listen,3)+size(all_ERP_talk,3)+1) = plot(EEG.times, grandaverage_ERP_listen(chani,:),'color', main_blue, 'LineWidth',2);
xlim([-250 750])
ylim([y_lim_upper y_lim_lower])
xlabel('Time [ms]')
title('ERPs for all Conditions')
ylabel('Amplitude [µV]')
p(1+size(all_ERP_listen,3)+size(all_ERP_talk,3)+2) = ...
    fill([50 150 150 50], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
legend(p([1 end-3 end-2 end-1]), 'talk - indiviual', 'listen - individual', 'talk - grand average', 'listen - grand average')
hold off

set(gca,'XTick',-250:50:750)

set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'grand_average_erp.png'))

%% End of processing

%display sanity check variables
ok_subj = cell2table(ok_subj, 'VariableNames',{'subj','subj_check_ID','time'})
writetable(ok_subj,fullfile(OUTPATH, '_erp_analysis_ok_subj.xlsx'))

check_done = 'OK'
save(fullfile(OUTPATH, '_erp_analysis_data.mat'), 'check_done')

save(fullfile(OUTPATH, '_erp_analysis_plot_data.mat'), 'grandaverage_ERP_talk', ...
    'grandaverage_ERP_listen', 'EEG', ...
    'all_ERP_listen', 'all_ERP_talk')

close(wb)
close all

