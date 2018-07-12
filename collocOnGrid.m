clc; clear;

%% Settings
grids = 1;
m = 20;
symmetric = 1;
kernel = 'gauss';
pde = 'square';

%% Setup
[rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde);
error = zeros(size(m));
amount_points = length(m);
[Xte, xlow, xup, ylow, yup] = collocation_points(w,34, grids);
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
    error(k) = max(abs(s_u - z));
    amount_points(k) = size(Xin,1);
    k = k + 1;
end
plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric)


%% Plot
figure
semilogy(amount_points, error)
xlabel('amount of collocation points')
ylabel('max. absolute error')
title('Wendland kernel')