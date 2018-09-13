clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'weighted'; % Choose weighted or standard collocation
grid = 1; % Collocation Points random or on grid
m = 20; % Amount and steps of Collocation points
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose the kernel
pde = 'square'; % Choose the PDE
error = 'abs'; % Calculate the error in the residuum or absolute


addpath(setting)
collocOnGrid(grid, m, symmetric, kernel, pde, error)
rmpath(setting)