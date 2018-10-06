function A_eval = evaluation_matrix(rbf, lap_rbf, gamma, Xin, Xte, w, symmetric)
if symmetric == 0
    A_eval = rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).',Xte(:,3),Xin(:,3).',Xte(:,4),Xin(:,4).').* w(Xte(:,1),Xte(:,2), Xte(:,3), Xte(:,4));
else
    A_eval = lap_rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).',Xte(:,3),Xin(:,3).',Xte(:,4),Xin(:,4).');
end
end