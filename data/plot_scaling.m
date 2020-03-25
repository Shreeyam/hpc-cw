clc;
clear;
close all;

x = [1 2 4];
y = [1823 5728 10246];

plot(x, y/1000, 'x--');
improvePlot;
grid;

xlabel('Cores [-]');
ylabel('Time Taken [s]');