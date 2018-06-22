clc; clear;

%% Settings
grid = 1;
m = 4:4:40;

%% Setup
[rbf, lap_rbf, f, w, realSol, realSolPlot] = allFunctions();
error = zeros(size(m));
amount_points = length(m);
Xte = collocation_points(w,10, grid);
grideval = collocation_points(w,100, grid);
z = realSolPlot(grideval(:,1), grideval(:,2));
k = 1;

%%
for i = m
    i
    Xin = collocation_points(w,i,grid);
    [gamma, alpha] = solvePDE(rbf, lap_rbf, w, Xin, Xte, f, realSol);
    A_eval = evaluation_matrix(rbf, gamma, Xin, grideval, w);
    s_u = A_eval*alpha;
    error(k) = max(abs(s_u - z));
    amount_points(k) = size(Xin,1);
    k = k + 1;
end



%% Plot
figure
semilogy(amount_points, error)
xlabel('amount of collocation points')
ylabel('max. absolute error')
title('Wendland kernel')