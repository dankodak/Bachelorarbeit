function Xbd = boundary_points(n, pde)
switch pde
    case "square"
        Xbd = [linspace(-1,1,n).',ones(n,1);
            linspace(-1,1,n).',-ones(n,1);
            ones(n,1), linspace(-1,1,n).';
            -ones(n,1),linspace(-1,1,n).'];
        Xbd = unique(Xbd, 'rows');
    case "circle"
        a = linspace(0,2*pi,n).';
        Xbd = [cos(a(1:n-1)), sin(a(1:n-1))];
    case "disc"
        a = linspace(0,2*pi,n).';
        Xbd = [cos(a(1:n-1)), sin(a(1:n-1))
            2*cos(a(1:n-1)), 2*sin(a(1:n-1))];
end

