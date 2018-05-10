clc; clear;
%Initialisierung der Punkte, Gitterweite, Epsilon und der Gleichungen
m = 15;
n = 100;
grid = 1;
rbf = @(eps,r) exp(-eps*r.^2);

% w = @(x) 1 - x(:,1).^2 - x(:,2).^2;
% f = @(x) exp(-(x(:,1).^2+x(:,2).^2)).*(-4+4*(x(:,1).^2+x(:,2).^2));
% % realSol = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
% realSol = @(x) 0;
% realSolPlot = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);

w = @(x) -x(:,1).^2 -x(:,2).^2 + x(:,1) + x(:,2) - sqrt(x(:,1).^4+x(:,2).^4-2*x(:,1).^3-2*x(:,2).^3+x(:,1).^2+x(:,2).^2);
f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSolPlot = @(x,y) sin(pi*x).*sin(pi*y);
realSolPlot = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% Bestimmung der Kollokations- und Testpunkte
[Xin, xlow, xup, ylow, yup] = collocation_points(w,m, grid);
Xte = test_points(xlow, xup, ylow, yup, m, w);

% Loese die PDE
[gamma, alpha] = solvePDE(rbf, w, Xin, Xte, f, realSol);


% Plot Kollokations- und Testpunkte
a = linspace(0,2*pi);
figure
axis equal
hold on
plot(Xin(:,1),Xin(:,2),'r+')% points inside
plot(Xte(:,1),Xte(:,2),'b*')
plot(cos(a),sin(a))
hold off


% Plot Loesung
[xx, yy] = ndgrid(linspace(xlow, xup, n),linspace(ylow, yup, n));
X = [xx(:), yy(:)];
bool = w(X) < 0;
A_eval = evaluation_matrix(rbf, gamma, Xin, X, w);
s_u = A_eval * alpha;
s_u(bool) = 0;
s_u = reshape(s_u,[n,n]);
z = realSolPlot(X);
z(bool) = 0;
z = reshape(z, [n,n]);

figure
subplot(2,2,1)
contour(xx,yy,s_u)
subplot(2,2,2)
contour(xx,yy,z)
subplot(2,2,3)
surf(xx,yy,s_u)
subplot(2,2,4)
surf(xx,yy,z)


maxerror = max(max(abs(s_u - z)))