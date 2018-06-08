function error = calculate_error(alpha, Xin, Xte, gamma, rbf, lap_rbf, f, w, realSol)
if realSol(0.5423, 0.823) == 0
    A = lap_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).');
    disc = A*alpha;
    real = f(Xte(:,1), Xte(:,2));
    error = max(abs(disc-real));
else
    Xte_eval = evaluation_matrix(rbf,gamma, Xin, Xte, w);
    app = Xte_eval*alpha;
    real = realSol(Xte(:,1), Xte(:,2));
    error = max(abs(app-real));
end