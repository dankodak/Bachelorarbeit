function error = calculate_error(alpha, Xin, Xte, gamma, rbf, lap_rbf, f, w, realSol)
if realSol([0.5,0.5]) == 0
    A = lap_rbf(gamma, Xte, Xin);
    disc = A*alpha;
    real = f(Xte);
    error = max(abs(disc-real));
else
    Xte_eval = evaluation_matrix(rbf,gamma, Xin, Xte, w);
    app = Xte_eval*alpha;
    real = realSol(Xte);
    error = max(abs(app-real));
end