clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
m = 15;
n = 100;
grid = 1;

[rbf, f, w, realSol, realSolPlot] = allFunctions();

% Bestimmung der Kollokations- und Testpunkte
[Xin, xlow, xup, ylow, yup] = collocation_points(w,m, grid);
Xte = test_points(xlow, xup, ylow, yup, m, w);

% Loese die PDE
[gamma, alpha] = solvePDE(rbf, w, Xin, Xte, f, realSol);


plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, n, realSolPlot)
% maxerror = max(max(abs(s_u - z)))