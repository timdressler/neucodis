% pp_plots.m
%
% Creates multiple plots for potential use for the poster/publication.
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

%% Plot: ERP and topography for talk and listen condition

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH, 'data\analysis_data\erp_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\plots\erp_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%variables to edit
CHAN = 'Cz';
ERP_FROM = 50; %only affects ERP topoplots, not the main ERP (see pp_erp_analysis.m)
ERP_TILL = 150; %only affects ERP topoplots, not the main ERP (see pp_erp_analysis.m)

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

%start eeglab
eeglab nogui
%load data
load(fullfile(INPATH, '_erp_analysis_plot_data.mat'))
%remove not needed channel
grandaverage_ERP_talk(20,:) = [];
grandaverage_ERP_listen(20,:) = [];
EEG = pop_select( EEG, 'rmchannel',{'IO'});

%preparations
%get channel ID
chani = find(strcmp({EEG.chanlocs.labels}, CHAN));
%get values for dynamic plot limits
y_lim_upper = min([min(grandaverage_ERP_talk(chani,:))  min(grandaverage_ERP_listen(chani,:)) ])-1;
y_lim_lower = max([max(grandaverage_ERP_talk(chani,:))  max(grandaverage_ERP_listen(chani,:)) ])+1;
%get ERP window
[~,win_start] = min(abs(EEG.times-ERP_FROM));
[~,win_end] = min(abs(EEG.times-ERP_TILL));

%create plot
%create main plot
close all
figure;
plot(EEG.times, grandaverage_ERP_talk(chani,:),'color', main_red, 'LineWidth',2);
hold on
plot(EEG.times, grandaverage_ERP_listen(chani,:),'color', main_blue, 'LineWidth',2);
xlim([-250 750])
ylim([y_lim_upper y_lim_lower])
% % xlim([-100 400])
xlabel('Time [ms]')
ylabel('Amplitude [µV]')
% % title('ERP for talk and listen condition')
fill([50 150 150 50], [y_lim_upper y_lim_upper y_lim_lower y_lim_lower], light_blue, 'FaceAlpha',0.1, 'EdgeColor','none');
l = legend('talk - grand average', 'listen - grand average');
fontsize(l,10,'points')
hold off

set(gca,'XTick',-250:50:750)

%add topoplot
%get limit values
cb_lim_lower = min([mean(grandaverage_ERP_listen(:,win_start:win_end),2), mean(grandaverage_ERP_talk(:,301:401),2)],[],'all');
cb_lim_upper = max([mean(grandaverage_ERP_listen(:,win_start:win_end),2), mean(grandaverage_ERP_talk(:,301:401),2)],[],'all');
%create topoplot
axes('Position',[.69 .16 .2 .2])
box on
topoplot(mean(grandaverage_ERP_talk(:,win_start:win_end),2), EEG.chanlocs, 'emarker2', {chani,'o','r',5,1})
title('N1 talk', 'Position', [0, 0.6, 0])
colormap("parula")
cb = colorbar;
title(cb, 'Amplitude [µV]')
clim([cb_lim_lower cb_lim_upper])

axes('Position',[.50 .16 .2 .2])
box on
topoplot(mean(grandaverage_ERP_listen(:,win_start:win_end),2), EEG.chanlocs, 'emarker2', {chani,'o','r',5,1})
title('N1 listen', 'Position', [0, 0.6, 0])
colormap("parula")
clim([cb_lim_lower cb_lim_upper])

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
%%saveas(gcf,fullfile(OUTPATH, 'pp_plot_erp_topo.png'))
exportgraphics(gcf,fullfile(OUTPATH, 'pp_plot_erp_topo.png'),'Resolution',1000)


%% Plot: dwPLI difference between talk and listen condition (Fronto-Temporal pairs)

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

%start eeglab
eeglab nogui

%load data
load(fullfile(INPATH, '_coherence_analysis_plot_data.mat'))

%get scale limits
for pairs = 1:4
    talk_GA_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1))), [], 'all'); 
    listen_GA_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1))), [], 'all'); 
    effect_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)) - squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1))), [], 'all'); 
end
scale_lim = max([talk_GA_scale, listen_GA_scale, effect_scale], [], 'all');

%create plots
close all

sp = 1;
figure;
for pairs = 1:2 %fronto-temporal 
    %extract needed values
    time = all_wpli(pairs).talk_GA.time; %time vector
    freq = all_wpli(pairs).talk_GA.freq; %frequency vector
    sig_clusters = squeeze(all_wpli(1).comparison.mask);  % Significance mask
    talk_GA_extracted = squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)); %grand average wPLI values (talk condition)
    listen_GA_extracted = squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1)); %grand average wPLI values (listen condition)
    effect = talk_GA_extracted - listen_GA_extracted; %wPLI difference between conditions
    stat = all_wpli(pairs).comparison; %cluster-based permutation test statistics

    %get significant clusters
    %positive clusters
    sig_pos_cluster = [];
    for clust = 1:length(stat.posclusters)
        if stat.posclusters(clust).prob < 0.05
            sig_pos_cluster(end+1) = clust;
        end
    end
    %positive clusters
    sig_neg_cluster = [];
    for clust = 1:length(stat.negclusters)
        if stat.negclusters(clust).prob < 0.05
            sig_neg_cluster(end+1) = clust;
        end
    end

    %get cluster locations
    %positive clusters
    pos_cluster_mat = squeeze(stat.posclusterslabelmat); %contains all clusters
    for tf = 1:numel(pos_cluster_mat) %removes all non-significant clusters
        if any(sig_pos_cluster == pos_cluster_mat(tf))
            %%pos_cluster_mat(tf) = 1;
        else
            pos_cluster_mat(tf) = 0;
        end
    end
    %negative clusters
    neg_cluster_mat = squeeze(stat.negclusterslabelmat); %contains all clusters
    for tf = 1:numel(neg_cluster_mat)
        if any(sig_neg_cluster == neg_cluster_mat(tf))
            %%neg_cluster_mat(tf) = 1;
        else
            neg_cluster_mat(tf) = 0;
        end
    end

    %plot wPLI listen condition
    subplot(2,3,sp)
    imagesc(time, freq, listen_GA_extracted);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Listen');
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI";
    caxis([-scale_lim scale_lim]);
    
    hold on
    text(0,0,'Test')
    hold off

    sp = sp+1;

    %plot wPLI talk condition
    subplot(2,3,sp)
    imagesc(time, freq, talk_GA_extracted);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Talk');
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI";
    caxis([-scale_lim scale_lim]);

    sp = sp+1;

    %plot wPLI difference between talk and listen condition with clusters highlighted (red = positive clusters, blue = negative clustera)
    subplot(2,3,sp)
    imagesc(time, freq, effect);
    caxis([-scale_lim scale_lim])
    axis xy;
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI difference";
    hold on;

    %overlay significant positive clusters with red contour
    if ~isempty(sig_pos_cluster)
        [C, h] = contour(time, freq, pos_cluster_mat, 1, 'r', 'LineWidth', 2);
    end

    %overlay significant negative clusters with blue contour
    if ~isempty(sig_neg_cluster)
        [C, h] = contour(time, freq, neg_cluster_mat, 1, 'b', 'LineWidth', 2);
    end

    hold off;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Difference: Talk - Listen');

    sp = sp+1;
end

%add row labels
labels = {'Fronto-Temporal L', 'Fronto-Temporal R'};
rows = 2;
for row = 1:rows
    ax = subplot(rows, 3, (row - 1) * 3 + 1); 
    %get the center y-position of the subplot 
    ax_pos = get(ax, 'Position');
    y_pos = ax_pos(2) + ax_pos(4) / 2 -0.1; 
    x_pos = 0.12; %shift label on x-axis
    annotation('textbox', [x_pos, y_pos - 0.05, 0.1, 0.1], ...
               'String', labels{row}, ...
               'FontSize', 14, ...
               'FontWeight', 'bold', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'Rotation', 90);
end

AddLetters2Plots('HShift', -0.02, 'VShift', -0.05, 'Location','NorthWest')


%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
%%saveas(gcf,fullfile(OUTPATH, 'pp_plot_dwpli_fronto_temporal.png'))
exportgraphics(gcf,fullfile(OUTPATH, 'pp_plot_dwpli_fronto_temporal.png'),'Resolution',1000)

% Plot: dwPLI difference between talk and listen condition (Fronto-Occipital pairs)

MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
INPATH = fullfile(MAINPATH, 'data\analysis_data\coherence_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
OUTPATH = fullfile(MAINPATH, 'data\plots\coherence_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)

%start eeglab
eeglab nogui

%load data
load(fullfile(INPATH, '_coherence_analysis_plot_data.mat'))

%get scale limits
for pairs = 1:4
    talk_GA_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1))), [], 'all'); 
    listen_GA_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1))), [], 'all'); 
    effect_scale(pairs) = max(abs(squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)) - squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1))), [], 'all'); 
end
scale_lim = max([talk_GA_scale, listen_GA_scale, effect_scale], [], 'all');

%create plots
close all

sp = 1;
figure;
for pairs = 3:4 %fronto-occipital 
    %extract needed values
    time = all_wpli(pairs).talk_GA.time; %time vector
    freq = all_wpli(pairs).talk_GA.freq; %frequency vector
    sig_clusters = squeeze(all_wpli(1).comparison.mask);  % Significance mask
    talk_GA_extracted = squeeze(mean(all_wpli(pairs).talk_GA.wpli_debiasedspctrm,1)); %grand average wPLI values (talk condition)
    listen_GA_extracted = squeeze(mean(all_wpli(pairs).listen_GA.wpli_debiasedspctrm,1)); %grand average wPLI values (listen condition)
    effect = talk_GA_extracted - listen_GA_extracted; %wPLI difference between conditions
    stat = all_wpli(pairs).comparison; %cluster-based permutation test statistics

    %get significant clusters
    %positive clusters
    sig_pos_cluster = [];
    for clust = 1:length(stat.posclusters)
        if stat.posclusters(clust).prob < 0.05
            sig_pos_cluster(end+1) = clust;
        end
    end
    %positive clusters
    sig_neg_cluster = [];
    for clust = 1:length(stat.negclusters)
        if stat.negclusters(clust).prob < 0.05
            sig_neg_cluster(end+1) = clust;
        end
    end

    %get cluster locations
    %positive clusters
    pos_cluster_mat = squeeze(stat.posclusterslabelmat); %contains all clusters
    for tf = 1:numel(pos_cluster_mat) %removes all non-significant clusters
        if any(sig_pos_cluster == pos_cluster_mat(tf))
            %%pos_cluster_mat(tf) = 1;
        else
            pos_cluster_mat(tf) = 0;
        end
    end
    %negative clusters
    neg_cluster_mat = squeeze(stat.negclusterslabelmat); %contains all clusters
    for tf = 1:numel(neg_cluster_mat)
        if any(sig_neg_cluster == neg_cluster_mat(tf))
            %%neg_cluster_mat(tf) = 1;
        else
            neg_cluster_mat(tf) = 0;
        end
    end

    %plot wPLI listen condition
    subplot(2,3,sp)
    imagesc(time, freq, listen_GA_extracted);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Listen');
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI";
    caxis([-scale_lim scale_lim]);
    
    hold on
    text(0,0,'Test')
    hold off

    sp = sp+1;

    %plot wPLI talk condition
    subplot(2,3,sp)
    imagesc(time, freq, talk_GA_extracted);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Talk');
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI";
    caxis([-scale_lim scale_lim]);

    sp = sp+1;

    %plot wPLI difference between talk and listen condition with clusters highlighted (red = positive clusters, blue = negative clustera)
    subplot(2,3,sp)
    imagesc(time, freq, effect);
    caxis([-scale_lim scale_lim])
    axis xy;
    colormap parula
    colorbar;
    cb=colorbar;
    cb.Title.String = "dwPLI difference";
    hold on;

    %overlay significant positive clusters with red contour
    if ~isempty(sig_pos_cluster)
        [C, h] = contour(time, freq, pos_cluster_mat, 1, 'r', 'LineWidth', 2);
    end

    %overlay significant negative clusters with blue contour
    if ~isempty(sig_neg_cluster)
        [C, h] = contour(time, freq, neg_cluster_mat, 1, 'b', 'LineWidth', 2);
    end

    hold off;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Difference: Talk - Listen');

    sp = sp+1;
end

%add row labels
labels = {'Fronto-Occipital L', 'Fronto-Occipital R'};
rows = 2;
for row = 1:rows
    ax = subplot(rows, 3, (row - 1) * 3 + 1); 
    %get the center y-position of the subplot 
    ax_pos = get(ax, 'Position');
    y_pos = ax_pos(2) + ax_pos(4) / 2 -0.1; 
    x_pos = 0.12; %shift label on x-axis
    annotation('textbox', [x_pos, y_pos - 0.05, 0.1, 0.1], ...
               'String', labels{row}, ...
               'FontSize', 14, ...
               'FontWeight', 'bold', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'Rotation', 90);
end

AddLetters2Plots('HShift', -0.02, 'VShift', -0.05, 'Location','NorthWest')


%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
%%saveas(gcf,fullfile(OUTPATH, 'pp_plot_dwpli_fronto_occipital.png'))
exportgraphics(gcf,fullfile(OUTPATH, 'pp_plot_dwpli_fronto_occipital.png'),'Resolution',1000)

%% Plot: TF-Plot with topographies
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

%load data
load(fullfile(INPATH, '_gamma_analysis_plot_data.mat'))

%start eeglab
eeglab nogui
%remove not needed channel
allersp_GRANDAVERAGE(:,:,20) = [];
allfreqs(:,:,20) = [];
alltimes(:,:,20) = [];
ersp_GRANDAVERAGE(:,20) = [];
EEG = pop_select( EEG, 'rmchannel',{'IO'});

%create plot
close all
figure;
%create topoplots
sub2 = subplot(2,3,4);
[~,time_idx1] = min(abs(alltimes(1,:,1)-(-130)));
[~,freq_idx1] = min(abs(allfreqs(1,:,1)-35));
topoplot(squeeze(allersp_GRANDAVERAGE(freq_idx1,time_idx1,:)),EEG.chanlocs, 'emarker2', {chani,'o','r',2,1})
title([num2str(-130) 'ms, ' num2str(35) 'Hz'],'Position', [0, 0.6, 0])

cb = colorbar;
colormap(parula);
clim([-1 1])
title(cb, 'Amplitude [dB]','Position', [0, 157, 0])

sub3 = subplot(2,3,5);
[~,time_idx2] = min(abs(alltimes(1,:,1)-(45)));
[~,freq_idx2] = min(abs(allfreqs(1,:,1)-35));
topoplot(squeeze(allersp_GRANDAVERAGE(freq_idx2,time_idx2,:)),EEG.chanlocs, 'emarker2', {chani,'o','r',2,1})
title([num2str(45) 'ms, ' num2str(35) 'Hz'],'Position', [0, 0.6, 0])
cb = colorbar;
colormap(parula);
clim([-1 1])
title(cb, 'Amplitude [dB]','Position', [0, 157, 0])

sub4 = subplot(2,3,6);
[~,time_idx3] = min(abs(alltimes(1,:,1)-(45)));
[~,freq_idx3] = min(abs(allfreqs(1,:,1)-45));
topoplot(squeeze(allersp_GRANDAVERAGE(freq_idx3,time_idx3,:)),EEG.chanlocs, 'emarker2', {chani,'o','r',2,1})
title([num2str(45) 'ms, ' num2str(45) 'Hz'],'Position', [0, 0.6, 0])
cb = colorbar;
colormap(parula);
clim([-1 1])
title(cb, 'Amplitude [dB]' ,'Position', [0, 157, 0])
%create TF plot
sub1 = subplot(2,3,1:3);
imagesc(ersp_times, ersp_freqs, ersp_GRANDAVERAGE)
ylabel('Hz')
xlabel('Time [ms]')
hold on
plot(ersp_times(time_idx1), ersp_freqs(freq_idx1), 'o', 'MarkerFaceColor','k')
text(ersp_times(time_idx1)+5, ersp_freqs(freq_idx1), 'B','color', 'k')
plot(ersp_times(time_idx2), ersp_freqs(freq_idx2), 'o', 'MarkerFaceColor','k')
text(ersp_times(time_idx2)+5, ersp_freqs(freq_idx2), 'C','color', 'k')
plot(ersp_times(time_idx3), ersp_freqs(freq_idx3), 'o', 'MarkerFaceColor','k')
text(ersp_times(time_idx3)+5, ersp_freqs(freq_idx3), 'D','color', 'k')
hold off
axis xy
xlim([-400 200])
ylim([15 60])
cb = colorbar;
colormap(parula);
clim([-1 1])
title(cb, 'Amplitude [dB]','Position', [0, 216, 0])

AddLetters2Plots({sub1, sub2, sub3, sub4}, {'A', 'B', 'C', 'D'},'HShift', -0.02, 'VShift', -0.05, 'Location','NorthWest')

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_tf_topo.png'))

%% Archive

%% Archive: Plot: TF-Plot with topographies and TF-topographies

% % MAINPATH = erase(SCRIPTPATH, 'neucodis\scripts');
% % INPATH = fullfile(MAINPATH, 'data\analysis_data\gamma_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
% % OUTPATH = fullfile(MAINPATH, 'data\plots\gamma_analysis'); %place 'data' folder in the same folder as the 'neucodis' folder %don't change names
% % FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
% % addpath(FUNPATH);
% % 
% % %sanity check
% % %check if folders exist
% % pp_check_folder_TD(MAINPATH, INPATH, OUTPATH)
% % %move current files to archive folder
% % % % pp_clean_up_folder_TD(OUTPATH)
% % 
% % %load data
% % load(fullfile(INPATH, '_gamma_analysis_plot_data.mat'))
% % 
% % %start eeglab
% % eeglab nogui
% % %remove not needed channel
% % allersp_GRANDAVERAGE(:,:,20) = [];
% % allfreqs(:,:,20) = [];
% % alltimes(:,:,20) = [];
% % ersp_GRANDAVERAGE(:,20) = [];
% % EEG = pop_select( EEG, 'rmchannel',{'IO'});
% % 
% % %create plot
% % % % close all
% % figure;
% % tftopo(allersp_GRANDAVERAGE,alltimes(:,:,1),allfreqs(:,:,1), ...
% %     'timefreqs', [-130 35; 45 35; 45 45], 'chanlocs', EEG.chanlocs, 'showchan', chani, 'limits', ...
% %     [-400 200 15 60 -1 1]);
% % % % sgtitle('Grand Average Topoplots');
% % colormap(parula);
% % 
% % %save plot
% % set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
% % % % saveas(gcf,fullfile(OUTPATH, 'pp_fig3_grand_average_tf_topo.png'))
% % 
% % %create plot
% % figure;
% % tftopo(allersp_GRANDAVERAGE,alltimes(:,:,1),allfreqs(:,:,1), 'chanlocs', ...
% %     EEG.chanlocs, 'showchan', chani, ...
% %     'plotscalponly', [-130 35]);
% % cb = colorbar;
% % colormap(parula);
% % clim([-1 1])
% % title(cb, 'Amplitude [dB]')
% % 
% % %save plot
% % set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
% % % % saveas(gcf,fullfile(OUTPATH, 'pp_fig4_grand_average_tf_topo.png'))
% % 
% % 
% % %end of processing
% % close all