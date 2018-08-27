function A_eval = evaluation_matrix(rbf, lap_rbf, gamma, Xin, Xbd, Xte, w, symmetric)
X = [Xin; Xbd];
if symmetric == 0
    A_eval = rbf(gamma, Xte(:,1), Xte(:,2), X(:,1).', X(:,2).');
else
    A_eval = [lap_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).'), rbf(gamma, Xte(:,1), Xte(:,2), Xbd(:,1).', Xbd(:,2).')];
end
end