%% Housekeeping

clc;
clear;
close all;

%% Load

vorticity = readmatrix('vorticity.csv');
stream = readmatrix('stream.csv');

%% Plot

disp(vorticity);
disp(stream);

figure;
sub