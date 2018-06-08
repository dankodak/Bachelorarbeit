function [Xin, xlow, xup, ylow, yup] = collocation_points(w,m,grid)

eps = 10^-10;

% Bounding Box
xlow = -1;
xup = 1;
ylow = -1;
yup = 1;

% check upper y bound
while true
    points = linspace(xlow, xup).';
    val = w(points, ones(size(points,1),1).*yup);
    if max(val) > 0 + eps
        yup = 2*yup;
        continue
    end
    break
end

% check lower y bound
while true
    points = linspace(xlow, xup).';
    val = w(points, ones(size(points,1),1).*ylow);
    if max(val) > 0 + eps
        ylow = 2*ylow;
        continue
    end
    break
end

% check upper x bound
while true
    points = linspace(ylow, yup).';
    val = w(ones(size(points,1),1).*xup, points);
    if max(val) > 0 + eps
        xup = 2*xup;
        continue
    end
    break
end

% check lower x bound
while true
    points = linspace(ylow, yup).';
    val = w(ones(size(points,1),1).*xlow, points);
    if max(val) > 0 + eps
        xlow = 2*xlow;
        continue
    end
    break
end


if grid == 1
    %Gitter erstellen
    [xx, yy] = ndgrid(linspace(xlow, xup, m),linspace(ylow,yup,m));
    X = [xx(:), yy(:)];

    %Bestimmen der Punkte im Inneren
    val = w(X(:,1), X(:,2));
    bool = val > 0 + eps;
    Xin = X(bool,:);
else
    Xin = zeros(m,2);
    k = 1;
    while Xin(m,:) == [0,0]
        pointx = xlow + (xup - xlow)*rand(1,1);
        pointy = ylow + (yup - ylow)*rand(1,1);
        point = [pointx pointy];
        val = w(point(1), point(2));
        if val > 0 + eps
            Xin(k,:) = point;
            k=k+1;
        end
    end
end

end