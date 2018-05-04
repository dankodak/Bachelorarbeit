function [Xin, Xbd, Nin, Nbd] = collocation_points(w,m)
eps = 10^-10

%Gitter erstellen
[xx, yy] = ndgrid(linspace(-1, 1, m));
X = [xx(:), yy(:)];

%Bestimmen der Punkte im Inneren
val = w(X(:,1), X(:,2));
bool = val > 0 + eps;
Xin = X(bool,:);

%Bestimmen der Punkte auf dem Rand
lin = linspace(0,2*pi,m);
Xbd = [sin(lin).',cos(lin).'];

[Nin,~] = size(Xin);
[Nbd,~] = size(Xbd);

% % Plot
% figure
% 
% axis equal
% 
% hold on
% plot(Xin(:,1),Xin(:,2),'r+') % points inside
% plot(Xbd(:,1),Xbd(:,2),'bo') % points outside
% hold off
end