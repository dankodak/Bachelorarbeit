function error = calculate_error(alpha, Xin, Xbd, Xte, gamma, rbf, lap_rbf, lap2_rbf, f, w, realSol, symmetric)
if realSol(0.5423, 0.823) == 0
    X = [Xin; Xbd];
    if symmetric == 0
        A = lap_rbf(gamma, Xte(:,1), Xte(:,2), X(:,1).', X(:,2).');
    else
        A = [lap2_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).'), lap_rbf(gamma, Xte(:,1), Xte(:,2), Xbd(:,1).', Xbd(:,2).')];
    end
    disc = A*alpha;
    real = f(Xte(:,1), Xte(:,2));
    error = max(abs(disc-real));
else
    Xte_eval = evaluation_matrix(rbf,lap_rbf, gamma, Xin, Xbd, Xte, w, symmetric);
    app = Xte_eval*alpha;
    real = realSol(Xte(:,1), Xte(:,2));
    error = max(abs(app-real));
end