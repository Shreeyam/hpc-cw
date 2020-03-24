clc;
clear;
close all;

load('fullA.mat');

b = [1:9]';
x = zeros(9, 1);

r0 = b - A * x;
p0 = r0;

for i = 1:3
ak = (dot(r0, r0))/dot(p0, A * p0);

x = x + ak .* p0;
rk1 = r0 - ak .* (A * p0);
eps = norm(rk1);
bk = dot(rk1, rk1) / dot(r0, r0);

r0 = rk1;

p0 = rk1 + bk * p0;
end
