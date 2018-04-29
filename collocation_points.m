function [Xin, Xbd, Nin, Nbd] = collocation_points(P,m)
xv = P(:,1);
yv = P(:,2);

%Gitter erstellen
[xx, yy] = ndgrid(linspace(0, 1, m));
X = [xx(:), yy(:)];
xq = X(:,1);
yq = X(:,2);

%Bestimmen der Punkte im Inneren
[in,on] = inpolygon(xq,yq,xv,yv);
in = xor(in,on);
Xin = X(in,:);

%Bestimmen der Punkte auf dem Rand
Xbd = zeros(1,2);
k = 1;
for i = 1:size(P)-1
    amount = floor(m*norm(P(i+1,:)-P(i,:)));
    dir = (P(i+1,:)-P(i,:))/amount;
    for j=1:amount
        Xbd(k,:) = P(i,:) + j*dir;
        k = k + 1;
    end
end

[Nin,~] = size(Xin);
[Nbd,~] = size(Xbd);

% Plot
% figure
% 
% plot(xv,yv) % polygon
% axis equal
% 
% hold on
% plot(xq(in),yq(in),'r+') % points inside
% plot(xq(~in),yq(~in),'bo') % points outside
% plot(Xbd(:,1),Xbd(:,2),'gd')
% hold off
end