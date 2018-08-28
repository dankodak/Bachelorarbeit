function greedy(n, symmetric, kernel, pde)
    % clc; clear;
    % warning off MATLAB:nearlySingularMatrix
    %% Settings
    % n = 30;
    % symmetric = 0;
    % kernel = 'wendland';
    % pde = 'square';

    %% Setup

    [rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
    error = zeros(1,n);
    gamma = zeros(1,n);

    % Bestimmung der Kollokations- und Testpunkte
    [Xin, xlow, xup, ylow, yup] = collocation_points(w,0, 1);
    Xbd = boundary_points(5);
    Xte = collocation_points(w,31, 1);
    grideval = collocation_points(w,100, 1);
    z = realSolPlot(grideval(:,1), grideval(:,2));

    %%
    for i = 1:n
        i
        % Loese die PDE
        [gamma(i), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xbd, Xte, f, realSol, symmetric);
        A_eval = evaluation_matrix(rbf, lap_rbf, gamma(i), Xin, Xbd, Xte, w, symmetric);
        s_u = A_eval*alpha;
        [error(i) , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(i), alpha, Xin, Xbd, grideval, z, symmetric, 'res');
        Xin(end+1,:) = grideval(index(1),:);
        grideval(index,:) = [];
        z(index) = [];
    end

    %% Plot
    [gamma(n+1), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xbd, Xte, f, realSol, symmetric);
    [error(n+1) , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(n+1), alpha, Xin, Xbd, grideval, z, symmetric, 'res');
    gamma = gamma(2:n+1);
    error = error(2:n+1);
    amount_points = (1:n);

    plot_sol(Xin, Xbd, Xte, xlow, xup, ylow, yup, w, f, rbf, lap_rbf, lap2_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
end