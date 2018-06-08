function [rbf, lap_rbf,  f, w, realSol, realSolPlot] = allFunctions()
syms a b c d gammas rbfs lap_rbfs ws prods fs realSols realSolPlots
% Kernel
%     function y = kernel(eps,r)
%         r = pdist2(Xte,Xin);
%         bool = r <= eps;
%         y = zeros(size(r));
%         y(bool) =((eps-r(bool)).^4 .*(4*r(bool) + eps))/20;
%     end
% rbf = @kernel;

rbfs(gammas,a,b,c,d) = exp(-gammas*((c-a)^2 + (d-b)^2));

% PDE auf Quadrat
% w = @(x) -x(:,1).^2 -x(:,2).^2 + 2 - sqrt(x(:,1).^4 + x(:,2).^4 - 2*x(:,1).^2 - 2*x(:,2).^2 +2);
% f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
% % realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSol = @(x) 0;
% realSolPlot = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% lap_rbf = @(eps,x,y) exp(-eps .* ((y(:,2).' - x(:,2)).^2 + (y(:,1).' - x(:,1)).^2)) .*(-4 + ( 4 .* (x(:,2) - x(:,2).^3).^2)./(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4).^(3./2) + ( 4 .* (x(:,1) - x(:,1).^3).^2)./(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4).^(3./2) + (2 - 6 .* x(:,2).^2)./ sqrt(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4)+ (2 - 6 .* x(:,1).^2)./ sqrt(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4) - 4 .* (-y(:,2).' + x(:,2)) .* eps .* (-2 .* x(:,2) - (2 .* x(:,2) .* (-1 + x(:,2).^2))./sqrt( 2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4))- 4 .* eps .* (-y(:,1).' + x(:,1)) .* (-2 .* x(:,1) - (2 .* x(:,1) .* (-1 + x(:,1).^2))./sqrt( 2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4)) + 4 .* eps .* (-2 + x(:,2).^2 + x(:,1).^2 + sqrt(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4)) - 4 .* (y(:,2).' - x(:,2)).^2 .* eps.^2 .* (-2 + x(:,2).^2 + x(:,1).^2 + sqrt(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4))- 4 .* eps.^2 .* (y(:,1).' - x(:,1)).^2 .* (-2 + x(:,2).^2 + x(:,1).^2 + sqrt(2 - 2 .* x(:,2).^2 + x(:,2).^4 - 2 .* x(:,1).^2 + x(:,1).^4)));



%PDE auf Kreis
ws(a,b) = 1 - a^2 - b^2;
fs(a,b) = exp(-a^2-b^2)*(-4+4*a^2+b^2);
realSols(a,b) = exp(-a^2-b^2) - 1/exp(1);
% realSols(a,b) = 0;
realSolPlots(a,b) = exp(-a^2-b^2) - 1/exp(1);
% realSolPlot = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);


rbf = matlabFunction(rbfs);
w = matlabFunction(ws);
f = matlabFunction(fs);
realSol = matlabFunction(realSols);
realSolPlot = matlabFunction(realSolPlots);
prods(gammas,a,b,c,d) = rbfs(gammas,a,b,c,d) * w(a,b);
lap_rbf = matlabFunction(diff(prods,a,2) + diff(prods,b,2));

end