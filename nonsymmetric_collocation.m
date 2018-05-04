clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
P = [0 0;0 1;1 1;1 0];
m = 15;
n = 100;
w = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);
f = @(x,y) exp(-(x.^2+y.^2)).*(-4+4*(x.^2+y.^2));
g = @(x,y) 0;
realSol = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);

%Bestimmung der Punkte, Ableitungen und der Kollokationsmatrix
[Xin, Xbd, Nin, Nbd] = collocation_points(w,m);
[rbf, lap_rbf] = RBFderivatives();
Xte = [Xin;Xbd];

[gamma, alpha] = solvePDE(rbf, lap_rbf, Xin, Xbd, Xte, Nin, Nbd, f, g, realSol);


[xx, yy] = ndgrid(linspace(-1, 1, n));
X = [xx(:), yy(:)];
A_eval = evaluation_matrix(rbf, gamma, Xin, Xbd, X);
s_u = A_eval * alpha;
s_u = reshape(s_u,[n,n]);
figure
subplot(2,1,1)
contour(xx,yy,s_u)
subplot(2,1,2)
contour(xx,yy,realSol(xx,yy))

% maxerror = max(max(abs(s_u - realSol(xx,yy))))
