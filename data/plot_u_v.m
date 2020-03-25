clc;
clear;
close all;

s_100 = csvread('streamfunc-161-100.csv');
s_400 = csvread('streamfunc-161-400.csv');
s_1000 = csvread('streamfunc-161-1000.csv');
s_3200 = csvread('streamfunc-161-3200.csv');

dx = 1/160;

u_100 = diff(s_100(81,:))/dx;
v_100 = -diff(s_100(:, 81))/dx;

u_400 = diff(s_400(81,:))/dx;
v_400 = -diff(s_400(:, 81))/dx;

u_1000 = diff(s_1000(81,:))/dx;
v_1000 = -diff(s_1000(:, 81))/dx;

u_3200 = diff(s_3200(81,:))/dx;
v_3200 = -diff(s_3200(:, 81))/dx;
%% u plots

figure;
hold on;
subplot(1, 2, 1);
plot(linspace(0, 1, 161), u_100(1:end-1)); hold on;
plot(linspace(0, 1, 161), u_400(1:end-1));
plot(linspace(0, 1, 161), u_1000(1:end-1));
plot(linspace(0, 1, 161), u_3200(1:end-1));
xlim([0 1]);
xlabel('L_y [m]')
ylabel('u [m/s]');
improvePlot;
xlim([0 1]);
legend('Re=100','Re=400', 'Re=1000', 'Re=3200');
grid;

%% v plots
subplot(1, 2, 2);
plot(linspace(0, 1, 160), v_100); hold on;
plot(linspace(0, 1, 160), v_400);
plot(linspace(0, 1, 160), v_1000);
plot(linspace(0, 1, 161), u_3200(1:end-1));

xlim([0 1]);
xlabel('L_x [m]')
ylabel('v [m/s]');
improvePlot;

legend('Re=100','Re=400', 'Re=1000', 'Re=3200');
grid;
