function [rbf, lap_rbf, lap2_rbf] = wendland(rbf1, lap_rbf1, lap2_rbf1)
    function y = kern(gamma,x1,x2,y1,y2)
        dist = pdist2([x1,x2],[y1.',y2.']);
        dist = dist.^2;
        bool = dist <= 1/gamma;
        y = rbf1(gamma,x1,x2,y1,y2) .* bool;
    end

    function y = laplace(gamma,x1,x2,y1,y2)
        dist = pdist2([x1,x2],[y1.',y2.']);
        dist = dist.^2;
        bool = dist <= 1/gamma;
        y = lap_rbf1(gamma,x1,x2,y1,y2) .* bool;
    end

    function y = laplace2(gamma,x1,x2,y1,y2)
        dist = pdist2([x1,x2],[y1.',y2.']);
        dist = dist.^2;
        bool = dist <= 1/gamma;        
        y = lap2_rbf1(gamma,x1,x2,y1,y2) .* bool;
    end
rbf = @kern;
lap_rbf = @laplace;
lap2_rbf = @laplace2;
end