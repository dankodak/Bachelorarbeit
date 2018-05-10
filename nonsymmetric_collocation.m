clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
m = 15;
n = 100;

w = @(x) 1 - x(:,1).^2 - x(:,2).^2;
f = @(x) exp(-(x(:,1).^2+x(:,2).^2)).*(-4+4*(x(:,1).^2+x(:,2).^2));
% realSol = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
realSol = @(x) 0
realSolPlot = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);

% w = @(x) -x(:,1).^2 -x(:,2).^2 + x(:,1) + x(:,2) - sqrt(x(:,1).^4+x(:,2).^4-2*x(:,1).^3-2*x(:,2).^3+x(:,1).^2+x(:,2).^2);
% f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSolPlot = @(x,y) sin(pi*x).*sin(pi*y);

%Bestimmung der Punkte, Ableitungen und der Kollokationsmatrix
[Xin, Nin] = collocation_points(w,m);
[rbf, lap_rbf] = RBFderivatives();
Xte = rand(m,2);

%Wende w auf die Stützstellen an
wonX = w(Xin);
wonX2 = w(Xte);
wonXeval = repmat(wonX2,[1,length(wonX)]);

[gamma, alpha] = solvePDE(rbf, w, Xin, Xte, f, realSol);


[xx, yy] = ndgrid(linspace(-1, 1, n));
X = [xx(:), yy(:)];
wonXplot = w(X);
wonXplot = repmat(wonXplot,[1,length(wonX)]);
A_eval = evaluation_matrix(rbf, gamma, Xin, X, w);
s_u = A_eval * alpha;
s_u = reshape(s_u,[n,n]);


figure
subplot(2,2,1)
contour(xx,yy,s_u)
subplot(2,2,2)
contour(xx,yy,realSolPlot(xx,yy))
subplot(2,2,3)
surf(xx,yy,s_u)
% axis([-1 1 -1 1 0 0.65])
subplot(2,2,4)
surf(xx,yy,realSolPlot(xx,yy))
% axis([-1 1 -1 1 0 0.65])

maxerror = max(max(abs(s_u - realSolPlot(xx,yy))))
