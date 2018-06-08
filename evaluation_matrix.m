function A_eval = evaluation_matrix(rbf, gamma, Xin, Xte, w)
%Bestimmen der Distance-Matrix
B = pdist2(Xte,Xin);

%Anwenden des Kernels
% A_eval = bsxfun(@times, rbf(gamma,B), w(Xte));
A_eval = bsxfun(@times,rbf(gamma, Xte(:,1), Xte(:,2), Xin(:,1).', Xin(:,2).'), w(Xte(:,1),Xte(:,2)));
end