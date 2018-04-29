function A_eval = evaluation_matrix(rbf, ep, Xin, Xbd, Xte)
%Bestimmen der Distance-Matrix
X = [Xin;Xbd];
B = pdist2(Xte,X);

%Anwenden des Kernels
A_eval = rbf(ep,B);
end