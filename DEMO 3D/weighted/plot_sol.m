function plot_sol(Xin, Xte, xlow, xup, ylow, yup, w, f, rbf, lap_rbf, lap2_rbf, gamma, alpha, realSolPlot, symmetric, amount_points, error)

%% Plot Fehler gegen Anzahl Kollokationspunkte
figure
semilogy(amount_points, error)
xlabel('amount of collocation points')
ylabel('max. error in derivative/absolute')
title('error plot')

%% Plot Bestes Gamma gegen Anzahl Kollokationspunkte
figure
semilogy(amount_points, gamma)
xlabel('amount of collocation points')
ylabel('gamma')
title('gamma plot')

end

