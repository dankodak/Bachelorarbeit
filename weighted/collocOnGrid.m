function collocOnGrid(grids, test_grid, m, test, symmetric, kernel, pde, calc_error)

    %% Setup
    [rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric);
    error = zeros(size(m));
    gamma = zeros(size(m));
    amount_points = length(m);
    [Xte, xlow, xup, ylow, yup] = collocation_points(w,test, test_grid);
    grideval = collocation_points(w,100, 1);
    z = realSolPlot(grideval(:,1), grideval(:,2));
    k = 1;

    %%
    for i = m
        i
        Xin = collocation_points(w,i,grids);
        [gamma(k), alpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric);
        A_eval = evaluation_matrix(rbf, lap_rbf, gamma(k), Xin, grideval, w, symmetric);
        s_u = A_eval*alpha;
        [error(k), ~] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma(k), alpha, Xin, grideval, z, symmetric, calc_error);
        amount_points(k) = size(Xin,1);
        k = k + 1;
    end
    plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, f, rbf, lap_rbf, lap2_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
end