function A_eval = evaluation_matrix(rbf, gamma, Xin, Xte, w)
%Bestimmen der Distance-Matrix
B = pdist2(Xte,Xin);

%Anwenden des Kernels
% A_eval = bsxfun(@times, rbf(gamma,B), w(Xte));
A_eval = bsxfun(@times,rbf(gamma, Xte, Xin), w(Xte));
end