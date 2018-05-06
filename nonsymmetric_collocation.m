clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
m = 20;
n = 100;
w = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);
f = @(x,y) exp(-(x.^2+y.^2)).*(-4+4*(x.^2+y.^2));
realSol = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);

%Bestimmung der Punkte, Ableitungen und der Kollokationsmatrix
[Xin, Nin] = collocation_points(w,m);
[rbf, lap_rbf] = RBFderivatives();
Xte = Xin;

%Wende w auf die Stützstellen an
wonX = w(Xin(:,1),Xin(:,2));
wonX2 = w(Xte(:,1),Xte(:,2));
wonXeval = repmat(wonX2,[1,length(wonX)]);

[gamma, alpha] = solvePDE(rbf, w, Xin, Xte, f);


[xx, yy] = ndgrid(linspace(-1, 1, n));
X = [xx(:), yy(:)];
wonXplot = w(X(:,1),X(:,2));
wonXplot = repmat(wonXplot,[1,length(wonX)]);
A_eval = evaluation_matrix(rbf, gamma, Xin, X, w);
s_u = A_eval * alpha;
s_u = reshape(s_u,[n,n]);


figure
subplot(2,2,1)
contour(xx,yy,s_u)
subplot(2,2,2)
contour(xx,yy,realSol(xx,yy))
subplot(2,2,3)
surf(xx,yy,s_u)
% axis([-1 1 -1 1 0 0.65])
subplot(2,2,4)
surf(xx,yy,realSol(xx,yy))
axis([-1 1 -1 1 0 0.65])

% maxerror = max(max(abs(s_u - realSol(xx,yy))))
