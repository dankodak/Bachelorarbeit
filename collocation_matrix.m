function A_Lambda = collocation_matrix(rbf, w, gamma, Xin)
h = 0.001;


%Erstellen der distance-Matrix
X = Xin;

xph = [X(:,1)+h,X(:,2)];
xmh = [X(:,1)-h,X(:,2)];
yph = [X(:,1),X(:,2)+h];
ymh = [X(:,1),X(:,2)-h];

wX = repmat(w(X),[1,length(X)]);
wxph = repmat(w(xph),[1,length(xph)]);
wxmh = repmat(w(xmh),[1,length(xph)]);
wyph = repmat(w(yph),[1,length(xph)]);
wymh = repmat(w(ymh),[1,length(xph)]);



B = squareform(pdist(X));
xphdist = pdist2(xph,X);
xmhdist = pdist2(xmh,X);
yphdist = pdist2(yph,X);
ymhdist = pdist2(ymh,X);



%Erstellen der Bloecke
Alap = (rbf(gamma, xmhdist).*wxmh + rbf(gamma, xphdist).*wxph + rbf(gamma, ymhdist).*wymh + rbf(gamma,yphdist).*wyph - 4*rbf(gamma , B).*wX)/h^2;


%Zusammensetzen
A_Lambda = Alap;
end