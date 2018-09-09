clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'weighted'; % Choose weighted or standard collocation
n = 30; % amount of iterations
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose your kernel
pde = 'square'; % Choose your PDE
error = "abs";


addpath(setting)
greedy(n, symmetric, kernel, pde, error);
rmpath(setting)