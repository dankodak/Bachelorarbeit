clc; clear;
grid = 1;
m = 30;
error = zeros(size(m));
k = 1;
for i = m
    i
    error(k) = nonsymmetric_collocation(i, grid);
    k = k + 1;
end
figure
semilogy(error)
xlabel('amount of collocation points')
ylabel('absolute error')
title('Wendland kernel')