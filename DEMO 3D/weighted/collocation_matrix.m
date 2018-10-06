function A_Lambda = collocation_matrix(lap_rbf, lap2_rbf, gamma, Xin, symmetric)
if symmetric == 0
    A_Lambda = lap_rbf(gamma, Xin(:,1), Xin(:,2), Xin(:,1).', Xin(:,2).',Xin(:,3), Xin(:,3).');
else
    A_Lambda = lap2_rbf(gamma, Xin(:,1), Xin(:,2), Xin(:,1).', Xin(:,2).',Xin(:,3), Xin(:,3).');
end
end