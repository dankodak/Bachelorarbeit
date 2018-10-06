clc; clear;
warning off MATLAB:nearlySingularMatrix

m = 1:5:21; % Amount and steps of Collocation points
m = m.^2;
test = 2500; % Amount of Test Points
symmetric = 0; % Symmetric or non symmetric collocation
error = 'abs'; % Calculate the error in the residuum or absolute

addpath('weighted')
collocOnGrid(0, 0, m, test, symmetric, 'gauss', 'square', error)
rmpath('weighted')