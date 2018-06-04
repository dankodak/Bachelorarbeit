function [rbf, lap_rbf,  f, w, realSol, realSolPlot] = allFunctions()
% Kernel
%     function y = kernel(eps,r)
%         r = pdist2(Xte,Xin);
%         bool = r <= eps;
%         y = zeros(size(r));
%         y(bool) =((eps-r(bool)).^4 .*(4*r(bool) + eps))/20;
%     end
% rbf = @kernel;
rbf = @(eps,x,y) exp(-eps*(bsxfun(@minus,x(:,1), y(:,1).').^2 + bsxfun(@minus,x(:,2), y(:,2).').^2));

% PDE auf Quadrat
w = @(x) -x(:,1).^2 -x(:,2).^2 + x(:,1) + x(:,2) - sqrt(x(:,1).^4+x(:,2).^4-2*x(:,1).^3-2*x(:,2).^3+x(:,1).^2+x(:,2).^2);
f = @(x) - 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
realSol = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
realSolPlot = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
    function y = laplace_gauss_quadrat(eps,x,y)
        e = exp(-eps*(bsxfun(@minus,x(:,1), y(:,1).').^2 + bsxfun(@minus,x(:,2), y(:,2).').^2));
        denom = 2 - 2 * x(:,2).^2 + x(:,2).^4 - 2 * x(:,1).^2 + x(:,1).^4;
        a1 = 4*(x(:,2) - x(:,2).^3).^2;
        a2 = 4*(x(:,1) - x(:,1).^3).^2;
        a3 = 2 - 6 * x(:,2).^2;
        a4 = 2 - 6 * x(:,1).^2;
        s1 = 4*(x(:,2)-y(:,2).') .* eps .*(-2*x(:,2) - (2*x(:,2).*(-1+x(:,2).^2)./sqrt(denom)));
        s2 = 4*(x(:,1)-y(:,1).') .* eps .*(-2*x(:,1) - (2*x(:,1).*(-1+x(:,1).^2)./sqrt(denom)));
        s3 = 2*eps*(-1+2*eps*(-x(:,2) + y(:,2).').^2).*(-2+x(:,2) + x(:,1) + sqrt(denom));
        s4 = 2*eps*(-1+2*eps*(-x(:,1) + y(:,1).').^2).*(-2+x(:,2) + x(:,1) + sqrt(denom));
        y = e .* (-4 + (a1 + a2)./((denom).^(3/2)) + (a3+a4)./sqrt(denom) - s1 - s2 - s3 - s4);
    end
lap_rbf = @laplace_gauss_quadrat;


%PDE auf Kreis
% w = @(x) 1 - x(:,1).^2 - x(:,2).^2;
% f = @(x) exp(-(x(:,1).^2+x(:,2).^2)).*(-4+4*(x(:,1).^2+x(:,2).^2));
% % realSol = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
% realSol = @(x) 0;
% realSolPlot = @(x) exp(-(x(:,1).^2+x(:,2).^2)) - 1/exp(1);
% lap_rbf = @(eps, x, y) -4*exp(-eps*(bsxfun(@minus,x(:,1), y(:,1).').^2 + bsxfun(@minus,x(:,2), y(:,2).').^2)).*(eps^2*(x(:,2).^2 + x(:,1).^2 - 1).*(bsxfun(@minus,x(:,1), y(:,1).').^2 + bsxfun(@minus,x(:,2), y(:,2).').^2) + eps*(bsxfun(@plus, 2 * x(:,1) * y(:,1).' + 2 * x(:,2) * y(:,2).' + 1, - 3 * x(:,2).^2 - 3 * x(:,1).^2))+1);

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