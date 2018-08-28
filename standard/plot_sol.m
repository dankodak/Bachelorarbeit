function plot_sol(Xin, Xbd, Xte, xlow, xup, ylow, yup, w, f, rbf, lap_rbf, lap2_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)
n = 100;


%% Plot Kollokations- und Testpunkte
% a = linspace(0,2*pi);
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

%% Plot Vergleich approximierte und analytische L�sung
figure
imagesc(abs(s_u - z))
colorbar
title('absolute error when compared to exact solution')

%% Plot Vergleich approximierte und analytische L�sung in der Ableitung
points = [Xin; Xbd];
if symmetric == 0
    A = lap_rbf(gamma(end), X(:,1), X(:,2), points(:,1).', points(:,2).');
else
    A = [lap2_rbf(gamma(end), X(:,1), X(:,2), Xin(:,1).', Xin(:,2).'), lap_rbf(gamma(end), X(:,1), X(:,2), Xbd(:,1).', Xbd(:,2).')];
end
disc = A*alpha;
disc(bool) = 0;
disc = reshape(disc,[n,n]);
real = f(X(:,1), X(:,2));
real(bool) = 0;
real = reshape(real,[n,n]);
figure
imagesc(abs(disc - real))
colorbar
title('error in derivative when compared to PDE')

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

