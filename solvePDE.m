function [bestgamma,retalpha] = solvePDE(rbf, lap_rbf, Xin, Xbd, Xte, Nin, Nbd, f, g, realSol)

gamma = -1.5:0.125:1.75;
gamma = 10.^gamma;

%Auswertung der DGL an den Teststellen
b = zeros(Nin+Nbd,1);
b(1:Nin) = f(Xte(1:Nin,1),Xte(1:Nin,2));
b(Nin+1:Nin+Nbd) = g(Xte(Nin+1:Nin+Nbd,1),Xte(Nin+1:Nin+Nbd,2));

%Echte Lösung
real = realSol(Xte(1:Nin+Nbd,1),Xte(1:Nin+Nbd,2));


error = zeros(size(gamma));
minerror = Inf;
bestgamma = gamma(1);

for i = 1:length(gamma)
    A_Lambda = collocation_matrix(rbf, lap_rbf, gamma(i), Xin, Xbd);

    %Bestimmung der Evaluationsmatrix
    A_eval = evaluation_matrix(rbf, gamma(i), Xin, Xbd, Xte);

    %Approximieren der DGL
    alpha = A_Lambda\b;
    s_u = A_eval * alpha;

    %Vergleich
    error(i) = max(abs((s_u - real)));
    if error(i) < minerror
        minerror = error(i);
        bestgamma = gamma(i);
        sol = s_u;
        retalpha = alpha;
    end
    
    
end
figure
semilogx(gamma, error);
% axis([-inf inf 0 1])
end