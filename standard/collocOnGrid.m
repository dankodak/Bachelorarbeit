clc; clear;
warning off MATLAB:nearlySingularMatrix
%% Settings
grids = 1;
m = 4:4:28;
symmetric = 0;
kernel = 'gauss';
pde = 'square';

%% Setup
[rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
error = zeros(size(m));
gamma = zeros(size(m));
amount_points = length(m);
[Xte, xlow, xup, ylow, yup] = collocation_points(w,31, grids);
grideval = collocation_points(w,100, grids);
z = realSolPlot(grideval(:,1), grideval(:,2));
k = 1;

%%
for i = m
    i
    Xin = collocation_points(w,i,grids);
    Xbd = boundary_points(i);
    [gamma(k), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xbd, Xte, f, realSol, symmetric);
    A_eval = evaluation_matrix(rbf, lap_rbf, gamma(k), Xin, Xbd, grideval, w, symmetric);
    s_u = A_eval*alpha;
    [error(k), ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(k), alpha, Xin, Xbd, grideval, z, symmetric, 'res');
    amount_points(k) = size(Xin,1);
    k = k + 1;
end
plot_sol(Xin, Xbd, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)