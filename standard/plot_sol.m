function plot_sol(Xin, Xbd, Xte, xlow, xup, ylow, yup, w, rbf, lap_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
n = 100;


%% Plot Kollokations- und Testpunkte
a = linspace(0,2*pi);
figure
axis equal
hold on
plot(Xin(:,1),Xin(:,2),'r+')% points inside
plot(Xbd(:,1),Xbd(:,2),'go')% points on border
plot(Xte(:,1),Xte(:,2),'b*')
% plot(cos(a),sin(a))
hold off


%% Plot Loesung
[xx, yy] = ndgrid(linspace(xlow, xup, n),linspace(ylow, yup, n));
X = [xx(:), yy(:)];
bool = w(X(:,1), X(:,2)) < 0;
A_eval = evaluation_matrix(rbf, lap_rbf, gamma(end), Xin, Xbd, X, w, symmetric);
s_u = A_eval * alpha;
s_u(bool) = 0;
s_u = reshape(s_u,[n,n]);
z = realSolPlot(X(:,1), X(:,2));
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

%% Plot Vergleich approximierte und analytische Lösung
figure
imagesc(abs(s_u - z))
colorbar

%% Plot Fehler gegen Anzahl Kollokationspunkte
figure
semilogy(amount_points, error)
xlabel('amount of collocation points')
ylabel('max. error in derivative/absolute')
title('Wendland/Gauss kernel')

%% Plot Bestes Gamma gegen Anzahl Kollokationspunkte
figure
semilogy(amount_points, gamma)
xlabel('amount of collocation points')
ylabel('gamma')
end

