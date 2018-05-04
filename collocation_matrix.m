function A_Lambda = collocation_matrix(rbf, lap_rbf, ep, Xin, Xbd)
h = 0.001;


%Erstellen der distance-Matrix
X = [Xin;Xbd];
[Nin,~] = size(Xin);
[Nbd,~] = size(Xbd);

xph = [X(:,1)+h,X(:,2)];
xmh = [X(:,1)-h,X(:,2)];
yph = [X(:,1),X(:,2)+h];
ymh = [X(:,1),X(:,2)-h];


B = squareform(pdist(X));
xphdist = pdist2(xph,X);
xmhdist = pdist2(xmh,X);
yphdist = pdist2(yph,X);
ymhdist = pdist2(ymh,X);

%Erstellen der Bloecke
% Alap = lap_rbf(ep , B(1:Nin,:));
Alap = (rbf(ep, xmhdist(1:Nin,:)) + rbf(ep, xphdist(1:Nin,:)) + rbf(ep, ymhdist(1:Nin,:)) + rbf(ep,yphdist(1:Nin,:)) - 4*rbf(ep , B(1:Nin,:)))/h^2;


Arbf = rbf(ep,B(Nin + 1 : Nin + Nbd,:));

%Zusammensetzen
A_Lambda = [Alap;Arbf];
end