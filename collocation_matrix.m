function A_Lambda = collocation_matrix(rbf, lap_rbf, ep, Xin, Xbd)
%Erstellen der distance-Matrix
X = [Xin;Xbd];
[Nin,~] = size(Xin);
[Nbd,~] = size(Xbd);
B = squareform(pdist(X));

%Erstellen der Bloecke
Alap = lap_rbf(ep , B(1:Nin,:));
Arbf = rbf(ep,B(Nin + 1 : Nin + Nbd,:));

%Zusammensetzen
A_Lambda = [Alap;Arbf];
end