clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'standard'; % Choose weighted or standard collocation
grid = 1; % Collocation Points random or on grid
test_grid = 1; % Test Points random or on grid
m = 10; % Amount and steps of Collocation points
test = 17; % Amount of Test Points
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose the kernel
pde = 'circle'; % Choose the PDE
error = 'abs'; % Calculate the error in the residuum or absolute


addpath(setting)
collocOnGrid(grid, test_grid, m, test, symmetric, kernel, pde, error)
rmpath(setting)