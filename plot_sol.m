function plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, rbf, gamma, alpha, realSolPlot)
n = 100;


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

end

