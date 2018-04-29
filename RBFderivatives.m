function [rbf, lap_rbf] = RBFderivatives()
rbf = @(eps,r) exp(-eps*r.^2);
lap_rbf = @(eps,r) exp(-eps*r.^2).*(-4*eps + 4*eps^2*r.^2);
end