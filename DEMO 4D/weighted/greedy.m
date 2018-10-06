function greedy(n, test_grid, test, symmetric, kernel, pde, calc_error)

    %% Setup

    [rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
    error = zeros(1,n);
    gamma = zeros(1,n);

    % Bestimmung der Kollokations- und Testpunkte
    [Xin, xlow, xup, ylow, yup] = collocation_points(w,0, 1);
    Xte = collocation_points(w,test, test_grid);
    grideval = collocation_points(w,100, 1);
    z = realSolPlot(grideval(:,1), grideval(:,2));

    %%
    for i = 1:n
        i
        % Loese die PDE
        [gamma(i), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
        % plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, realSolPlot)
        A_eval = evaluation_matrix(rbf, lap_rbf, gamma(i), Xin, grideval, w, symmetric);
        [error(i) , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(i), alpha, Xin, grideval, z, symmetric, 'res');
        if calc_error == "abs"
            [error(i) , ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(i), alpha, Xin, grideval, z, symmetric, 'abs');
        end
        Xin(end+1,:) = grideval(index(1),:);
        grideval(index,:) = [];
        z(index) = [];
    end

    %% Plot
    [gamma(n+1), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
    if calc_error == "abs"
        [error(n+1) , ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(n+1), alpha, Xin, grideval, z, symmetric, 'abs');
    else
        [error(n+1) , ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(n+1), alpha, Xin, grideval, z, symmetric, 'res');
    end
    gamma = gamma(2:n+1);
    error = error(2:n+1);
    amount_points = (1:n);

    plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, f, rbf, lap_rbf, lap2_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
end