clc;
clear;
close all;

dx = 0.25;
dy = 0.25;
dx2 = dx^2;
dy2 = dy^2;

Ny = 5;
Nx = 5;

A = diag(ones(1, Ny * Nx) * (2/(dx2) + 2/(dy2))) + ... 
    diag(ones(1, Ny * Nx - 1) * -1/dy2, 1) + diag(ones(1, Ny * Nx - 1) * -1/dy2, -1) + ...
    diag(ones(1, (Ny-1) * Nx) * -1/dy2, Ny) + diag(ones(1, (Ny-1) * Nx) * -1/dy2, -Ny);


v = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 -0.0128 -0.0128205 -0.0128205 0; 0 -8 -8 -8 0];
v_vec = reshape(v, [Nx * Ny, 1]);

psi = inv(A) * v_vec;

psi_matrix = reshape(psi, [Ny, Nx]);
disp(psi_matrix);