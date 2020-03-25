%% Housekeeping

clc;
clear;
% close all;

%% Load

vorticity = csvread('dumpmatrix.csv');
% stream = readmatrix('stream.csv');

%% Plot

% disp(vorticity);
% disp(stream);

figure;
contourf(vorticity);
improvePlot;
title("C++ Sim Vorticity");
% sub