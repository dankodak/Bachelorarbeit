clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'weighted'; % Choose weighted or standard collocation
grids = 1; % Collocation Points random or on grid
m = 4:4:28; % Amount and steps of Collocation points
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose the kernel
pde = 'newone'; % Choose the PDE
error = 'res'; % Calculate the error in the residuum or absolute


addpath(setting)
collocOnGrid(grids, m, symmetric, kernel, pde, error)
rmpath(setting)