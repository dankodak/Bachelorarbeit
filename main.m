clc; clear;
grid = 1;
m = 50;
error = zeros(size(m));
for i = m
    error(i) = nonsymmetric_collocation(i, grid);
end
figure
semilogy(error)
xlabel('amount of collocation points')
ylabel('absolute error')
title('Wendland kernel')