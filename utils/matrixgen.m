clc;
clear;
close all;

Lx = 1;
Ly = 1;

Ny = 100;
Nx = 100;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dx2 = dx^2;
dy2 = dy^2;

U = 1;
Re = 100;
dt = (Re * dx * dy)/4;

subpdiags = ones(1, (Ny-2) * (Nx-2) - 1) * (-1/dy2);
subpdiags(Ny-2:Ny-2:end) = 0;

A = diag(ones(1, (Ny-2) * (Nx-2)) * (2/(dx2) + 2/(dy2))) + ... 
    diag(subpdiags, 1) + diag(subpdiags, -1) + ...
    diag(ones(1, (Ny-3) * (Nx-2)) * -1/dy2, Ny-2) + diag(ones(1, (Ny-3) * (Nx-2)) * -1/dy2, -(Ny-2));

A = sparse(A);

v = zeros(Ny, Nx);
s = zeros(Ny, Nx);
% v = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 -0.0128 -0.0128205 -0.0128205 0; 0 -8 -8 -8 0];
% v_vec = reshape(v, [Nx * Ny, 1]);
% 
% psi = inv(A) * v_vec;
% 
% s = reshape(psi, [Ny, Nx]);
% s0 = s;
% disp(s);

%% Implementing algorithm in MATLAB

for i = 1:1
% Update vorticity boundary conditions

% Top
v(Ny, :) = (s(Ny, :) - s(Ny-1, :)) * (2/dy2) - (2 * U)/dy;

% Bottom
v(1, :) = (s(1, :) - s(2, :)) * (2/dy2);

% Left
v(:, 1) = (s(:, 1) - s(:, 2)) * (2/dx2);

% Right
v(:, Nx) = (s(:, Nx) - s(:, Nx-1)) * (2/dx2);

% Calculate interior vorticity
for i = 2:Nx-1
    for j = 2:Ny-1
        v(j, i) = -(((s(j, i+1) - 2 * s(j, i) + s(j, i-1))/(dx2)) + ((s(j+1, i) - 2 * s(j, i) + s(j-1, i))/(dy2)));
    end
end

% Calculate new vorticity

for i = 2:Nx-1
    for j = 2:Ny-1
        term1 = -1 * (((s(j+1, i) - s(j-1, i))/(2 * dy))*((v(j, i+1) - v(j, i-1))/(2 * dx)));
        term2 = (((s(j, i+1) - s(j, i-1))/(2 * dx))*((v(j+1, i) - v(j-1, i))/(2 * dy)));
        term3 = (1/Re) * (((v(j, i+1) - 2 * v(j, i) + v(j, i-1))/(dx2))+((v(j+1, i) - 2 * v(j, i) + v(j-1, i))/(dy2)));
        
        v(j, i) = v(j, i) + dt * (term1 + term2 + term3);
    end
end

v_vec = reshape(v(2:end-1, 2:end-1), [(Nx-2) * (Ny-2), 1]);
psi_vec = A\v_vec;
s(2:end-1, 2:end-1) = reshape(psi_vec, [Ny-2, Nx-2]);
end

figure;
imagesc(v);
improvePlot;
title("MATLAB Sim Vorticity");

%%
figure;
vorticity = csvread('dumpmatrix.csv');
imagesc(abs(v - vorticity(:, 1:end-2)));
title('Deltas');
colorbar;
