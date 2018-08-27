function [error , index] = greedy_error(rbf, lap_rbf, lap2_rbf, w, f, gamma, alpha, Xin, Xbd, grideval, z, symmetric, setting)
switch setting
    case 'abs'
        A_eval = evaluation_matrix(rbf, lap_rbf, gamma, Xin, Xbd, grideval, w, symmetric);
        s_u = A_eval*alpha;
        error = max(abs(s_u - z));
        index = find(abs(s_u-z)==error);
    case 'res'
        X = [Xin; Xbd];
        if symmetric == 0
            A = lap_rbf(gamma, grideval(:,1), grideval(:,2), X(:,1).', X(:,2).');
        else
            A = [lap2_rbf(gamma, grideval(:,1), grideval(:,2), Xin(:,1).', Xin(:,2).'), lap_rbf(gamma, grideval(:,1), grideval(:,2), Xbd(:,1).', Xbd(:,2).')];
        end
        disc = A*alpha;
        real = f(grideval(:,1), grideval(:,2));
        error = max(abs(disc-real));
        index = find(abs(disc-real)==error);
end