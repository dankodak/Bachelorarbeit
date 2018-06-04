function A_Lambda = collocation_matrix(rbf, lap_rbf, w, gamma, Xin)
%Erstellen der distance-Matrix
B = squareform(pdist(Xin));

%Erstellen der Bloecke
A_Lambda = lap_rbf(gamma, B, Xin, Xin);
end