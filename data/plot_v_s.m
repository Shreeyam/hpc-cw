clc;
clear;
close all;

v = csvread('vorticity-161-100.csv');
s = csvread('streamfunc-161-100.csv');

figure;
subplot(1, 2, 1);
contourf(v);
colorbar;
improvePlot;
xlabel('L_x [m]');
ylabel('L_y [m]');
title('Vorticity \omega');
axis equal;
subplot(1, 2, 2);
contourf(s);
xlabel('L_x [m]');
ylabel('L_y [m]');
title('Streamfunction \Psi');
axis equal;
colorbar;
improvePlot;