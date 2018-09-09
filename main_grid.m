clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'standard'; % Choose weighted or standard collocation
grids = 1; % Collocation Points random or on grid
m = 25; % Amount and steps of Collocation points
symmetric = 1; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose the kernel
pde = 'square'; % Choose the PDE
error = 'abs'; % Calculate the error in the residuum or absolute


addpath(setting)
collocOnGrid(grids, m, symmetric, kernel, pde, error)
rmpath(setting)