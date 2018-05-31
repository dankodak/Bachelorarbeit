function A_Lambda = collocation_matrix(rbf, w, gamma, Xin)
h = 0.001;

%Erstellen der distance-Matrix

xph = [Xin(:,1)+h,Xin(:,2)];
xmh = [Xin(:,1)-h,Xin(:,2)];
yph = [Xin(:,1),Xin(:,2)+h];
ymh = [Xin(:,1),Xin(:,2)-h];

B = squareform(pdist(Xin));
xphdist = pdist2(xph,Xin);
xmhdist = pdist2(xmh,Xin);
yphdist = pdist2(yph,Xin);
ymhdist = pdist2(ymh,Xin);

%Erstellen der Bloecke
A_Lambda = (bsxfun(@times, rbf(gamma, xmhdist), w(xmh)) + bsxfun(@times, rbf(gamma, xphdist), w(xph)) + bsxfun(@times, rbf(gamma, ymhdist), w(ymh)) + bsxfun(@times, rbf(gamma, yphdist), w(yph)) - 4*bsxfun(@times, rbf(gamma , B), w(Xin)))/h^2;

end