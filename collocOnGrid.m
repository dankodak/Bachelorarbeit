clc; clear;

%% Settings
grids = 1;
m = 15;
symmetric = 0;
kernel = 'gauss';
pde = 'square';

%% Setup
[rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
error = zeros(size(m));
amount_points = length(m);
[Xte, xlow, xup, ylow, yup] = collocation_points(w,31, grids);
grideval = collocation_points(w,100, grids);
z = realSolPlot(grideval(:,1), grideval(:,2));
k = 1;

%%
for i = m
    i
    Xin = collocation_points(w,i,grids);
    [gamma, alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
    A_eval = evaluation_matrix(rbf, lap_rbf, gamma, Xin, grideval, w, symmetric);
    s_u = A_eval*alpha;
    [error(k), ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma, alpha, Xin, grideval, z, symmetric, 'res');
    amount_points(k) = size(Xin,1);
    k = k + 1;
end
plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric)


%% Plot
figure
semilogy(amount_points, error)
xlabel('amount of collocation points')
ylabel('max. error')
title('Wendland kernel')