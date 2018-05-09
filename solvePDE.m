function [bestgamma,retalpha] = solvePDE(rbf, w, Xin, Xte, f, realSol)

gamma = -1.5:0.125:3;
gamma = 10.^gamma;

%Auswertung der DGL an den Teststellen
b = f(Xin);

error = zeros(size(gamma));
condition = zeros(size(gamma));
minerror = Inf;
bestgamma = gamma(1);

for i = 1:length(gamma)
    A_Lambda = collocation_matrix(rbf, w, gamma(i), Xin);

    %Approximieren der DGL
    alpha = A_Lambda\b;
    
    %Vergleich
    error(i) = calculate_error(alpha, Xin, Xte, gamma(i), rbf, f, w, realSol);
    condition(i) = cond(A_Lambda);
    if error(i) < minerror
        minerror = error(i);
        bestgamma = gamma(i);
        retalpha = alpha;
    end
    
    
end
figure
subplot(1,2,1)
semilogx(gamma, error);
subplot(1,2,2)
semilogx(gamma, condition);
% axis([-inf inf 0 1])
end