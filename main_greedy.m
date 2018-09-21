clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'weighted'; % Choose weighted or standard collocation
n = 100; % amount of iterations
test_grid = 1; % Test Points random or on grid
test = 17; % Amount of Test Points
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose your kernel
pde = 'circle'; % Choose your PDE
error = "abs";


addpath(setting)
greedy(n, test_grid, test, symmetric, kernel, pde, error);
rmpath(setting)