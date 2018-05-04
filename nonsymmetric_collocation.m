clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
P = [0 0;0 1;1 1;1 0];
m = 12;
n = 50;
f = @(x,y) - 2*pi^2*sin(pi*x).*sin(pi*y);
g = @(x,y) 0;
realSol = @(x,y) sin(pi*x).*sin(pi*y);

%Bestimmung der Punkte, Ableitungen und der Kollokationsmatrix
[Xin, Xbd, Nin, Nbd] = collocation_points(P,m);
[rbf, lap_rbf] = RBFderivatives();
Xte = [Xin;Xbd];

[gamma, alpha] = solvePDE(rbf, lap_rbf, Xin, Xbd, Xte, Nin, Nbd, f, g, realSol);


[xx, yy] = ndgrid(linspace(0, 1, n));
X = [xx(:), yy(:)];
A_eval2 = evaluation_matrix(rbf, gamma, Xin, Xbd, X);
s_u2 = A_eval2 * alpha;
s_u2 = reshape(s_u2,[n,n]);
figure
subplot(2,1,1)
contour(xx,yy,s_u2)
subplot(2,1,2)
contour(xx,yy,realSol(xx,yy))

maxerror = max(max(abs(s_u2 - realSol(xx,yy))))
