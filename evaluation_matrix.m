function A_eval = evaluation_matrix(rbf, gamma, Xin, Xbd, Xte)
%Bestimmen der Distance-Matrix
X = [Xin;Xbd];
B = pdist2(Xte,X);

%Anwenden des Kernels
A_eval = rbf(gamma,B);
end