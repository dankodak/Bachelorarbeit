function A_eval = evaluation_matrix(rbf, gamma, Xin, Xte, w)
%Bestimmen der Distance-Matrix
X = Xin;
B = pdist2(Xte,X);

wonXte = repmat(w(Xte),[1,size(Xin,1)]);
%Anwenden des Kernels
A_eval = rbf(gamma,B);
A_eval = A_eval .* wonXte;
end