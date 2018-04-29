clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
P = [0 0;0 1;1 1;1 0];
m = 50;
f = @(x,y) - 2*pi^2*sin(pi*x).*sin(pi*y);
g = @(x,y) 0;
realSol = @(x,y) sin(pi*x).*sin(pi*y);

%Bestimmung der Punkte, Ableitungen und der Kollokationsmatrix
[Xin, Xbd, Nin, Nbd] = collocation_points(P,m);
[rbf, lap_rbf] = RBFderivatives();
Xte = [Xin;Xbd];

sol = solvePDE(rbf, lap_rbf, Xin, Xbd, Xte, Nin, Nbd, f, g, realSol);


% [xx, yy] = ndgrid(linspace(0, 1, m));
% X = [xx(:), yy(:)];
% A_eval2 = evaluation_matrix(rbf, ep, Xin, Xbd, X);
% s_u2 = A_eval2 * alpha;
% s_u2 = reshape(s_u2,[m,m]);
% subplot(2,1,1)
% surf(xx,yy,s_u2)
% subplot(2,1,2)
% surf(xx,yy,g(xx,yy))


% [xx,yy] = meshgrid(Xte(:,1),Xte(:,2));
% a = size(xx)
% contourf(xx, yy, reshape(sol, [a,a]))


%pcolor(Xte(:,1),Xte(:,2),sol)
