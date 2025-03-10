% pp_methods_illustration.m
%
% Creates plots illustrating the methods used.
% Using simulated data.
%
% Tim Dressler, 10.03.2025

clear all
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
OUTPATH = fullfile(MAINPATH, 'data\plots\methods_illustration'); 
FUNPATH = fullfile(MAINPATH, 'neucodis\functions\');
addpath(FUNPATH);

%sanity check
%check if folders exist
pp_check_folder_TD(MAINPATH, OUTPATH)
%move current files to archive folder
pp_clean_up_folder_TD(OUTPATH)

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

%set fontsize
font_size = 28;


%% Plot: The concept of phase (1)

%define parameters
fs = 1000; %sampling frequency (Hz)
T = 1; %duration (s)
t = linspace(0, T, fs*T); %time vector

%generate sine wave
s = sin(2*pi*t);

%define three time points for phase extraction
t_points = [0.1, 0.45, 0.75]; %example time points
colors = {main_blue, main_green, light_blue}; %colors for each point

%compute phase angles
phase_angles = 2 * pi * t_points;

%plot sine wave
figure;
plot(t, s, 'k', 'LineWidth', 1.5); % Plot full sine wave
hold on;
for i = 1:length(t_points)
    [~, idx] = min(abs(t - t_points(i))); % Find closest index
    xline(t(idx), '--', 'Color', colors{i}, 'LineWidth',1.5)
end
xlabel('Time (s)');
ylabel('Amplitude');
grid off;

%add legend
lgd = legend({'', [num2str(t_points(1)) 's'], [num2str(t_points(2)) 's'], [num2str(t_points(3)) 's']})
%%fontsize(lgd,26,'points')
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phase1.png'))

%% Plot: The concept of phase (2)

%create a polar axis for the phase plot
figure;
polar_ax = polaraxes;

%create the unit circle using polarplot
theta = linspace(0, 2*pi, 100);
polarplot(polar_ax, theta, ones(1,100), 'k--'); % Unit circle

%plot angular lines (spokes) and radial ticks
polar_ax.ThetaTick = 0:45:360; 
polar_ax.RLim = [0 1]; 
polar_ax.RTick = [1]; 

%plot lines for each time point on the polar plot
hold on;
for i = 1:length(t_points)
    %phase angle for the current time point
    phase = phase_angles(i);
    %plot line at the phase angle (radius is 1 because it's a unit circle)
    polarplot(polar_ax, [0 phase], [0 1], 'Color', colors{i}, 'LineWidth', 2); % Line at each phase
end

%add legend
lgd = legend({'', [num2str(t_points(1)) 's'], [num2str(t_points(2)) 's'], [num2str(t_points(3)) 's']});
%%fontsize(lgd,26,'points')
%add title 
title('Phase Angles at Difference Timepoints');
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phase2.png'))

%% Plot: Phase differences (1.1)

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

%define parameters
fs = 1000; %sampling frequency (Hz)
T = 1; %duration (s)
t = linspace(0, T, fs*T); %time vector

%generate first sine wave (without phase shift)
s1 = sin(2*pi*t);

%generate second sine wave with a phase shift
phi_shift = -pi/8; %phase shift 
s2 = sin(2*pi*t + phi_shift);

%define a time point for phase extraction
t_point = 0.25; % Example time point
colors = {main_blue, main_green}; % Colors for each wave

%compute phase angles at the specified time point
[~, idx] = min(abs(t - t_point)); % Find closest index for t_point
phase1 = angle(exp(1i * 2 * pi * t(idx))); % Phase of the first sine wave
phase2 = angle(exp(1i * (2 * pi * t(idx) + phi_shift))); % Phase of the second sine wave

%plot both sine waves
figure;
plot(t, s1,'Color', main_blue, 'LineWidth', 1.5); % Plot first sine wave
hold on;
plot(t, s2, 'Color', main_yellow, 'LineWidth', 1.5); % Plot second sine wave
xline(t(idx), '--', 'Color', 'k', 'LineWidth',1.5)
xlabel('Time (s)');
ylabel('Amplitude');
grid off;

%add legend
lgd = legend({'F8', 'T8'});
%%fontsize(lgd,26,'points')
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phasediff1_1.png'))

%% Plot: Phase differences (2.1)

%create a polar axis for the phase plot
figure;
polar_ax = polaraxes;

%create the unit circle using polarplot
theta = linspace(0, 2*pi, 100);
polarplot(polar_ax, theta, ones(1,100), 'k--'); 

%plot angular lines (spokes) and radial ticks
polar_ax.ThetaTick = 0:45:360; %set angular ticks (you can adjust the step)
polar_ax.RLim = [0 1]; %limit the radius to 1 (for the unit circle)
polar_ax.RTick = [1]; %show radial tick only at 1 (unit circle)

%plot lines for each sine wave phase at the given time point
hold on;
polarplot(polar_ax, [0 phase1], [0 1], 'Color', main_blue, 'LineWidth', 2); % Phase of first wave
polarplot(polar_ax, [0 phase2], [0 1], 'Color', main_yellow, 'LineWidth', 2); % Phase of second wave
polarplot(polar_ax, [0 phase1-phase2], [0 1], 'Color', 'r', 'LineWidth', 2); % Phase of second wave

%add legend
lgd = legend({'', 'F8', 'T8', 'Difference'});
%%fontsize(lgd,26,'points')
%add title 
title('Phase Angles of Both Signals and Phase Difference');
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phasediff2_1.png'))

%% Plot: Phase differences (1.2)

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

%define parameters
fs = 1000; %sampling frequency (Hz)
T = 1; %duration (s)
t = linspace(0, T, fs*T); %time vector

%generate first sine wave (without phase shift)
s1 = sin(2*pi*t);

%generate second sine wave with a phase shift
phi_shift = pi/8; %phase shift %changes compared to above
s2 = sin(2*pi*t + phi_shift);

%define a time point for phase extraction
t_point = 0.25; % Example time point
colors = {main_blue, main_green}; % Colors for each wave

%compute phase angles at the specified time point
[~, idx] = min(abs(t - t_point)); % Find closest index for t_point
phase1 = angle(exp(1i * 2 * pi * t(idx))); % Phase of the first sine wave
phase2 = angle(exp(1i * (2 * pi * t(idx) + phi_shift))); % Phase of the second sine wave

%plot both sine waves
figure;
plot(t, s1,'Color', main_blue, 'LineWidth', 1.5); % Plot first sine wave
hold on;
plot(t, s2, 'Color', main_yellow, 'LineWidth', 1.5); % Plot second sine wave
xline(t(idx), '--', 'Color', 'k', 'LineWidth',1.5)
xlabel('Time (s)');
ylabel('Amplitude');
grid off;

%add legend
lgd = legend({'F8', 'T8'});
%%fontsize(lgd,26,'points')
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phasediff1_2.png'))

%% Plot: Phase differences (2.2)

%create a polar axis for the phase plot
figure;
polar_ax = polaraxes;

%create the unit circle using polarplot
theta = linspace(0, 2*pi, 100);
polarplot(polar_ax, theta, ones(1,100), 'k--'); 

%plot angular lines (spokes) and radial ticks
polar_ax.ThetaTick = 0:45:360; %set angular ticks (you can adjust the step)
polar_ax.RLim = [0 1]; %limit the radius to 1 (for the unit circle)
polar_ax.RTick = [1]; %show radial tick only at 1 (unit circle)

%plot lines for each sine wave phase at the given time point
hold on;
polarplot(polar_ax, [0 phase1], [0 1], 'Color', main_blue, 'LineWidth', 2); % Phase of first wave
polarplot(polar_ax, [0 phase2], [0 1], 'Color', main_yellow, 'LineWidth', 2); % Phase of second wave
polarplot(polar_ax, [0 phase1-phase2], [0 1], 'Color', 'r', 'LineWidth', 2); % Phase of second wave

%add legend
lgd = legend({'', 'F8', 'T8', 'Difference'});
%%fontsize(lgd,26,'points')
%add title 
title('Phase Angles of Both Signals and Phase Difference');
fontsize(font_size,"points")

%save plot
set(gcf, 'Position', get(0, 'Screensize')-[0 0 300 150]);
saveas(gcf,fullfile(OUTPATH, 'pp_plot_phasediff2_2.png'))



