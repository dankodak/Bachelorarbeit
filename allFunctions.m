function [rbf, lap_rbf,  f, w, realSol, realSolPlot] = allFunctions()
syms a b c d gammas rbfs lap_rbfs ws prods fs realSols realSolPlots
%% Kernel
%     function y = kernel(eps,r)
%         r = pdist2(Xte,Xin);
%         bool = r <= eps;
%         y = zeros(size(r));
%         y(bool) =((eps-r(bool)).^4 .*(4*r(bool) + eps))/20;
%     end
% rbf = @kernel;

rbfs(gammas,a,b,c,d) = exp(-gammas*((c-a)^2 + (d-b)^2));

%% PDE auf Quadrat
ws(a,b)= -a^2 - b^2 + 2 -sqrt(a^4 - 2*a^2 + b^4 - 2*b^2 +2);
fs(a,b) = - 2*pi^2*sin(pi*a)*sin(pi*b);
realSols(a,b) = sin(pi*a)*sin(pi*b);
% realSols(a,b) = 0;
realSolPlots(a,b) = sin(pi*a)*sin(pi*b);



%% PDE auf Kreis
% ws(a,b) = 1 - a^2 - b^2;
% fs(a,b) = exp(-a^2-b^2)*(-4+4*(a^2+b^2));
% % realSols(a,b) = exp(-a^2-b^2) - 1/exp(1);
% realSols(a,b) = 0;
% realSolPlots(a,b) = exp(-a^2-b^2) - 1/exp(1);

%%
rbf = matlabFunction(rbfs);
w = matlabFunction(ws);
f = matlabFunction(fs);
realSol = matlabFunction(realSols);
realSolPlot = matlabFunction(realSolPlots);
prods(gammas,a,b,c,d) = rbfs(gammas,a,b,c,d) * w(a,b);
lap_rbf = matlabFunction(diff(prods,a,2) + diff(prods,b,2));

end