function error = calculate_error(alpha, Xin, Xte, gamma, rbf, lap_rbf, lap2_rbf, f, w, realSol, symmetric)
if realSol(0.5423, 0.823) == 0
    if symmetric == 0
        A = lap_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).',Xte(:,3), Xin(:,3).');
    else
        A = lap2_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).',Xte(:,3), Xin(:,3).');
    end
    disc = A*alpha;
    real = f(Xte(:,1), Xte(:,2), Xte(:,3));
    error = max(abs(disc-real));
else
    Xte_eval = evaluation_matrix(rbf,lap_rbf, gamma, Xin, Xte, w, symmetric);
    app = Xte_eval*alpha;
    real = realSol(Xte(:,1), Xte(:,2),Xte(:,3), Xte(:,4));
    error = max(abs(app-real));
end