function A_Lambda = collocation_matrix(rbf, lap_rbf, lap2_rbf, gamma, Xin, Xbd, symmetric)
X = [Xin; Xbd];
if symmetric == 0
    A_Lambda = [lap_rbf(gamma, Xin(:,1), Xin(:,2), X(:,1).', X(:,2).');
        rbf(gamma, Xbd(:,1), Xbd(:,2), X(:,1).', X(:,2).')];
else
    A_Lambda1 = lap2_rbf(gamma, Xin(:,1), Xin(:,2), Xin(:,1).', Xin(:,2).');
    A_Lambda2 = lap_rbf(gamma, Xin(:,1), Xin(:,2), Xbd(:,1).', Xbd(:,2).');
    A_Lambda3 = rbf(gamma, Xbd(:,1), Xbd(:,2), Xbd(:,1).', Xbd(:,2).');
    A_Lambda = [A_Lambda1, A_Lambda2; A_Lambda2.', A_Lambda3];
end
end