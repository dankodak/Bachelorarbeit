function [maxerror] = nonsymmetric_collocation(m, grid)

    [rbf, lap_rbf, f, w, realSol, realSolPlot] = allFunctions();

    % Bestimmung der Kollokations- und Testpunkte
    [Xin, xlow, xup, ylow, yup] = collocation_points(w,m, grid);
    Xte = test_points(xlow, xup, ylow, yup, m, w);

    % Loese die PDE
    [gamma, alpha] = solvePDE(rbf, lap_rbf, w, Xin, Xte, f, realSol);


    plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, realSolPlot)
    z = realSolPlot(Xin(:,1), Xin(:,2));
    A_eval = evaluation_matrix(rbf, gamma, Xin, Xin, w);
    s_u = A_eval*alpha;
    maxerror = max(max(abs(s_u - z)));
end