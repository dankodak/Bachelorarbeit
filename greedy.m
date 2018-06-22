clc; clear;

%% Settings
grid = 1;
m = 4;
n = 100;

%% Setup

[rbf, lap_rbf, f, w, realSol, realSolPlot] = allFunctions();
error = zeros(size(n));

% Bestimmung der Kollokations- und Testpunkte
[Xin, xlow, xup, ylow, yup] = collocation_points(w,m, grid);
Xte = collocation_points(w,10, grid);
grideval = collocation_points(w,100, grid);
z = realSolPlot(grideval(:,1), grideval(:,2));

%%
for i = 1:n
    i
    % Loese die PDE
    [gamma, alpha] = solvePDE(rbf, lap_rbf, w, Xin, Xte, f, realSol);
    % plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, realSolPlot)
    A_eval = evaluation_matrix(rbf, gamma, Xin, grideval, w);
    s_u = A_eval*alpha;
    error(i) = max(abs(s_u - z));
    index = find(abs(s_u-z)==error(i));
    Xin(end+1,:) = grideval(index,:);
    grideval(index,:) = [];
    z(index) = [];
end


%% Plot
figure
semilogy((1:n)+m,error)
xlabel('amount of collocation points')
ylabel('max. absolute error')
title('Wendland kernel')


figure
axis equal
hold on
plot(Xin(:,1),Xin(:,2),'r+')% points inside
plot(Xte(:,1),Xte(:,2),'b*')
hold off