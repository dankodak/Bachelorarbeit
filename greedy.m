clc; clear;
warning off MATLAB:nearlySingularMatrix
%% Settings
grid = 1;
m = 2;
n = 70;
symmetric = 0;
kernel = 'wendland';
pde = 'square';

%% Setup

[rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde);
error = zeros(size(n));
gamma = zeros(size(n));

% Bestimmung der Kollokations- und Testpunkte
[Xin, xlow, xup, ylow, yup] = collocation_points(w,m, grid);
Xte = collocation_points(w,34, grid);
grideval = collocation_points(w,100, grid);
z = realSolPlot(grideval(:,1), grideval(:,2));

%%
for i = 1:n
%     i
    % Loese die PDE
    [gamma(i), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
    % plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, realSolPlot)
    A_eval = evaluation_matrix(rbf, lap_rbf, gamma(i), Xin, grideval, w, symmetric);
    s_u = A_eval*alpha;
    [error(i) , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(i), alpha, Xin, grideval, z, symmetric, 'res');
    Xin(end+1,:) = grideval(index(1),:);
    grideval(index,:) = [];
    z(index) = [];
end


%% Plot
figure
semilogy((1:n)+m,error)
xlabel('amount of collocation points')
ylabel('max. absolute error')
title(kernel)

figure
semilogy((1:n)+m,gamma)
xlabel('amount of collocation points')
ylabel('gamma')
title(kernel)

% figure
% axis equal
% hold on
% plot(Xin(:,1),Xin(:,2),'r+')% points inside
% plot(Xte(:,1),Xte(:,2),'b*')
% hold off

[gamma, alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric)