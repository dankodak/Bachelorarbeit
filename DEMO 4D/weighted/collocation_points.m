function [Xin, xlow, xup, ylow, yup] = collocation_points(w,m,grid)

eps = 10^-10;

% Bounding Box
xlow = -1;
xup = 1;
ylow = -1;
yup = 1;
zlow = -1;
zup = 1;
klow = -1;
kup = 1;

% check upper y bound

if grid == 1
    %Gitter erstellen
    [xx, yy] = ndgrid(linspace(xlow, xup, m),linspace(ylow,yup,m));
    X = [xx(:), yy(:)];

    %Bestimmen der Punkte im Inneren
    val = w(X(:,1), X(:,2));
    bool = val > 0 + eps;
    Xin = X(bool,:);
else
    Xin = zeros(m,4);
    k = 1;
    while Xin(m,:) == [0,0,0,0]
        pointx = xlow + (xup - xlow)*rand(1,1);
        pointy = ylow + (yup - ylow)*rand(1,1);
        pointz = zlow + (zup - zlow)*rand(1,1);
        pointk = klow + (kup - klow)*rand(1,1);
        point = [pointx pointy pointz pointk];
        val = w(point(1), point(2), point(3), point(4));
        if val > 0 + eps
            Xin(k,:) = point;
            k=k+1;
        end
    end
end

end