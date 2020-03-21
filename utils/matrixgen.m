clc;
clear;
close all;

dx = 1;
dy = 1;
dx2 = dx^2;
dy2 = dy^2;

Ny = 3;
Nx = 3;

A = diag(ones(1, Ny * Nx) * (2/(dx2) + 2/(dy2))) + ... 
    diag(ones(1, Ny * Nx - 1) * -1/dy2, 1) + diag(ones(1, Ny * Nx - 1) * -1/dy2, -1) + ...
    diag(ones(1, (Ny-1) * Nx) * -1/dy2, Ny) + diag(ones(1, (Ny-1) * Nx) * -1/dy2, -Ny);

