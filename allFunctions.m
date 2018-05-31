function [rbf,  f, w, realSol, realSolPlot] = allFunctions()
% Kernel
    function y = kernel(eps,r)
        bool = r <= eps;
        y = zeros(size(r));
        y(bool) =((eps-r(bool)).^4 .*(4*r(bool) + eps))/20;
    end
rbf = @kernel;
% rbf = @(eps,r) exp(-eps*r.^2);

% w = @(x) -x(:,1).^2 -x(:,2).^2 + x(:,1) + x(:,2) - sqrt(x(:,1).^4+x(:,2).^4-2*x(:,1).^3-2*x(:,2).^3+x(:,1).^2+x(:,2).^2);
% f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSolPlot = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));


%PDE auf Kreis
w = @(x) 1 - x(:,1).^2 - x(:,2).^2;
f = @(x) exp(-(x(:,1).^2+x(:,2).^2)).*(-4+4*(x(:,1).^2+x(:,2).^2));
realSol = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
% realSol = @(x) 0;
realSolPlot = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);




end



% w = @(x) 1 - x(:,1).^2 - x(:,2).^2;
% f = @(x) exp(-(x(:,1).^2+x(:,2).^2)).*(-4+4*(x(:,1).^2+x(:,2).^2));
% % realSol = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
% realSol = @(x) 0;
% realSolPlot = @(x,y) exp(-(x.^2+y.^2)) - 1/exp(1);

% w = @(x) -x(:,1).^2 -x(:,2).^2 + x(:,1) + x(:,2) - sqrt(x(:,1).^4+x(:,2).^4-2*x(:,1).^3-2*x(:,2).^3+x(:,1).^2+x(:,2).^2);
% f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
% realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
% % realSolPlot = @(x,y) sin(pi*x).*sin(pi*y);
% realSolPlot = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));