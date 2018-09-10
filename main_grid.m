clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'standard'; % Choose weighted or standard collocation
grid = 1; % Collocation Points random or on grid
m = 4:1:45; % Amount and steps of Collocation points
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose the kernel
pde = 'square'; % Choose the PDE
error = 'res'; % Calculate the error in the residuum or absolute


addpath(setting)
collocOnGrid(grid, m, symmetric, kernel, pde, error)
rmpath(setting)