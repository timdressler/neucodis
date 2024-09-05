%% pp_fastfft_TD.m
clear
close all
clc

fs = 500;         
t1 = 1;           
t2 = 1;        
f1 = 5;          
f2 = 10;           

t_1 = 0:1/fs:t1-1/fs;   
t_2 = t1:1/fs:t1+t2-1/fs; 

y1 = sin(2*pi*f1*t_1);  %5Hz
y2 = sin(2*pi*f2*t_2);  %10Hz 

y = [y1 y2];


pp_fastfft_TD(y, 500, 1, 'wins', 1, 'overlap', 0.99, 'fig', 1, 'ylim_min', 0, 'ylim_max', 15)

