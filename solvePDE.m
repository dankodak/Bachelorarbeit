function [bestgamma,retalpha] = solvePDE(rbf, lap_rbf, lap2_rbf, w, Xin, Xte, f, realSol, symmetric)
step = 0.125;
gammapot = -1:step:2.5;
gamma = 10.^gammapot;

%Auswertung der DGL an den Teststellen
b = f(Xin(:,1), Xin(:,2));

error = zeros(size(gamma));
condition = zeros(size(gamma));
minerror = Inf;
bestgamma = gamma(1);

for i = 1:length(gamma)
    % Kollokationsmatrix erstellen
%     if symmetric == 0
%         A_Lambda = lap_rbf(gamma(i), Xin(:,1), Xin(:,2), Xin(:,1).', Xin(:,2).');
%     else
%         A_Lambda = lap2_rbf(gamma(i), Xin(:,1), Xin(:,2), Xin(:,1).', Xin(:,2).');
%     end
    A_Lambda = collocation_matrix(lap_rbf, lap2_rbf, gamma(i), Xin, symmetric);
    % LGS lösen
    alpha = A_Lambda\b;
    
    %Vergleich
    error(i) = calculate_error(alpha, Xin, Xte, gamma(i), rbf, lap_rbf, lap2_rbf, f, w, realSol, symmetric);
%     condition(i) = cond(A_Lambda);
    if error(i) < minerror
        minerror = error(i);
        bestgamma = gamma(i);
        retalpha = alpha;
    end
end

while bestgamma == gamma(1)
    gammapot = [(gammapot(1) - step) , gammapot];
    gamma = [10^gammapot(1) gamma];
    A_Lambda = collocation_matrix(lap_rbf, lap2_rbf, gamma(1), Xin, symmetric);

    %Approximieren der DGL
    alpha = A_Lambda\b;
    error = [calculate_error(alpha, Xin, Xte, gamma(1), rbf, lap_rbf, lap2_rbf, f, w, realSol, symmetric) error];
    condition = [cond(A_Lambda) condition];
    if error(1) < minerror
        minerror = error(1);
        bestgamma = gamma(1);
        retalpha = alpha;
    end
end

while bestgamma == gamma(end)
    gammapot = [gammapot , (gammapot(end) + step)];
    gamma = [gamma 10^gammapot(end)];
    A_Lambda = collocation_matrix(lap_rbf, lap2_rbf, gamma(end), Xin, symmetric);

    %Approximieren der DGL
    alpha = A_Lambda\b;
    
    error = [error calculate_error(alpha, Xin, Xte, gamma(1), rbf, lap_rbf, lap2_rbf, f, w, realSol, symmetric)];
    condition = [condition cond(A_Lambda)];
    if error(end) < minerror
        minerror = error(end);
        bestgamma = gamma(end);
        retalpha = alpha;
    end
end


minerror


% figure
% subplot(1,2,1)
% loglog(gamma, error);
% xlabel('gamma')
% ylabel('error')
% subplot(1,2,2)
% loglog(gamma, condition);
% xlabel('gamma')
% ylabel('condition')
end