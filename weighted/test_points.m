function Xte = test_points(xlow, xup, ylow, yup, n, w)
Xte = zeros(n,2);
k = 1;
while Xte(n,:) == [0,0]
    pointx = xlow + (xup - xlow)*rand(1,1);
    pointy = ylow + (yup - ylow)*rand(1,1);
    point = [pointx pointy];
    val = w(point(1), point(2));
    if val > 0 + eps
        Xte(k,:) = point;
        k=k+1;
    end
end

end