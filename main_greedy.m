clc; clear;
warning off MATLAB:nearlySingularMatrix


setting = 'weighted'; % Choose weighted or standard collocation
n = 100; % amount of iterations
symmetric = 0; % Symmetric or non symmetric collocation
kernel = 'gauss'; % Choose your kernel
pde = 'newone'; % Choose your PDE


addpath(setting)
greedy(n, symmetric, kernel, pde);
rmpath(setting)