function [Xin, Nin] = collocation_points(w,m)
eps = 10^-10;

%Gitter erstellen
[xx, yy] = ndgrid(linspace(-1, 1, m));
X = [xx(:), yy(:)];

%Bestimmen der Punkte im Inneren
val = w(X(:,1), X(:,2));
bool = val > 0 + eps;
Xin = X(bool,:);


[Nin,~] = size(Xin);

% Plot
% figure
% 
% axis equal
% 
% plot(Xin(:,1),Xin(:,2),'r+') % points inside

end