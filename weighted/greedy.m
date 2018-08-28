function greedy(n, symmetric, kernel, pde)
    % clc; clear;
    % warning off MATLAB:nearlySingularMatrix
    % %% Settings
    % grid = 1;
    % n = 10;
    % symmetric = 0;
    % kernel = 'gauss';
    % pde = 'square';

    %% Setup

    [rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
    error = zeros(1,n);
    gamma = zeros(1,n);

    % Bestimmung der Kollokations- und Testpunkte
    [Xin, xlow, xup, ylow, yup] = collocation_points(w,0, 1);
    Xte = collocation_points(w,31, 1);
    grideval = collocation_points(w,100, 1);
    z = realSolPlot(grideval(:,1), grideval(:,2));

    %%
    for i = 1:n
        i
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
    [gamma(n+1), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
    [error(n+1) , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(n+1), alpha, Xin, grideval, z, symmetric, 'res');
    gamma = gamma(2:n+1);
    error = error(2:n+1);
    amount_points = (1:n);

    plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
end