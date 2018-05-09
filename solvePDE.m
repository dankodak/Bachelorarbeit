function [bestgamma,retalpha] = solvePDE(rbf, w, Xin, Xte, f)

gamma = -1.5:0.125:3;
gamma = 10.^gamma;

%Auswertung der DGL an den Teststellen
% b = f(Xin(:,1),Xin(:,2));
b = f(Xin);

error = zeros(size(gamma));
minerror = Inf;
bestgamma = gamma(1);

for i = 1:length(gamma)
    A_Lambda = collocation_matrix(rbf, w, gamma(i), Xin);

    %Approximieren der DGL
    alpha = A_Lambda\b;
    
    %Vergleich
    error(i) = calculate_error(alpha, Xin, Xte, gamma(i), rbf, f, w);
    if error(i) < minerror
        minerror = error(i);
        bestgamma = gamma(i);
        retalpha = alpha;
    end
    
    
end
figure
semilogx(gamma, error);
% axis([-inf inf 0 1])
end