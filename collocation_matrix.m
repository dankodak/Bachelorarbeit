function A_Lambda = collocation_matrix(rbf, w, gamma, Xin)
h = 0.001;


%Erstellen der distance-Matrix
X = Xin;

xph = [X(:,1)+h,X(:,2)];
xmh = [X(:,1)-h,X(:,2)];
yph = [X(:,1),X(:,2)+h];
ymh = [X(:,1),X(:,2)-h];

wX = repmat(w(X),[1,length(X)]);
wxph = repmat(w(xph),[1,size(xph,1)]);
wxmh = repmat(w(xmh),[1,size(xmh,1)]);
wyph = repmat(w(yph),[1,size(yph,1)]);
wymh = repmat(w(ymh),[1,size(ymh,1)]);



B = squareform(pdist(X));
xphdist = pdist2(xph,X);
xmhdist = pdist2(xmh,X);
yphdist = pdist2(yph,X);
ymhdist = pdist2(ymh,X);

%Erstellen der Bloecke
Alap = (rbf(gamma, xmhdist).*wxmh + rbf(gamma, xphdist).*wxph + rbf(gamma, ymhdist).*wymh + rbf(gamma,yphdist).*wyph - 4*rbf(gamma , B).*wX)/h^2;
% Alap = lap_rbf(ep , B(1:Nin,:));

%Zusammensetzen
A_Lambda = Alap;
end