function Xbd = boundary_points(n)
Xbd = [linspace(-1,1,n).',ones(n,1);
    linspace(-1,1,n).',-ones(n,1);
    ones(n,1), linspace(-1,1,n).';
    -ones(n,1),linspace(-1,1,n).'];
Xbd = unique(Xbd, 'rows');
end

