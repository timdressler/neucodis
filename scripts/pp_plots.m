% pp_plots.m
%
% Creates the plot shown on the poster/in the publication.
% Using analysis data.
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

%% Plot 1, ERP for talk and listen condition

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH, 'data\analysis_data\erp_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\plots\erp_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
CHAN = 'Cz';

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%load data
load(fullfile(INPATH, '_erp_analysis_plot_data.mat'))

%preparations
%get channel ID
chani = find(strcmp({EEG.chanlocs.labels}, CHAN));
%get values for dynamic plot limits
y_lim_upper = min([min(grandaverage_ERP_talk(chani,:))  min(grandaverage_ERP_listen(chani,:)) ])-1;
y_lim_lower = max([max(grandaverage_ERP_talk(chani,:))  max(grandaverage_ERP_listen(chani,:)) ])+1;

%create plot
close all
figure;
plot(EEG.times, grandaverage_ERP_talk(chani,:),'color', main_red, 'LineWidth',2);
hold on
plot(EEG.times, grandaverage_ERP_listen(chani,:),'color', main_blue, 'LineWidth',2);
xlim([-250 750])
ylim([y_lim_upper y_lim_lower])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
title('ERP for talk and listen condition')
fill([50 150 150 50], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
legend('talk - grand average', 'listen - grand average')
hold off

set(gca,'XTick',-250:50:750)

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_fig1_grand_average_erp.png'))

%% Plot 2, dwPLI difference between talk and listen condition

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH, 'data\analysis_data\coherence_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\plots\coherence_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%load data



%% Plot 3, TF-Plot with topographies

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH, 'data\analysis_data\gamma_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\plots\gamma_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%start eeglab 
eeglab nogui

%load data
load(fullfile(INPATH, '_gamma_analysis_plot_data.mat'))

%create plot
figure;
tftopo(allersp_GRANDAVERAGE,alltimes(:,:,1),allfreqs(:,:,1), ...
    'timefreqs', [-130 45; 45 45; 45 40], 'chanlocs', EEG.chanlocs, 'showchan', chani);
sgtitle('Grand Average Topoplots');
colormap(parula);

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_fig3_grand_average_tf_topo.png'))




