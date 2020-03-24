%% Housekeeping

clc;
clear;
% close all;

%% Load

vorticity = csvread('161x161-1000steps.csv');
% stream = readmatrix('stream.csv');

%% Plot

% disp(vorticity);
% disp(stream);

figure;
contourf(vorticity);
improvePlot;
title("C++ Sim Vorticity");
% sub