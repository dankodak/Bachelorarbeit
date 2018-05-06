function error = calculate_error(alpha, Xin, Xte, gamma, rbf, f, w)
h = 0.001;


xph = [Xte(:,1)+h,Xte(:,2)];
xmh = [Xte(:,1)-h,Xte(:,2)];
yph = [Xte(:,1),Xte(:,2)+h];
ymh = [Xte(:,1),Xte(:,2)-h];

Xte_eval = evaluation_matrix(rbf,gamma, Xin, Xte, w);
xph_eval = evaluation_matrix(rbf,gamma, Xin, xph, w);
xmh_eval = evaluation_matrix(rbf,gamma, Xin, xmh, w);
yph_eval = evaluation_matrix(rbf,gamma, Xin, yph, w);
ymh_eval = evaluation_matrix(rbf,gamma, Xin, ymh, w);

Xte_sol = Xte_eval * alpha;
xph_sol = xph_eval * alpha;
xmh_sol = xmh_eval * alpha;
yph_sol = yph_eval * alpha;
ymh_sol = ymh_eval * alpha;

disc = (xph_sol + xmh_sol + yph_sol + ymh_sol - 4*Xte_sol)/h^2;

real = f(Xte(:,1), Xte(:,2));

error = max(abs(disc-real));
end