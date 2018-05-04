function error = calculate_error(alpha, Xin, Xbd, Xte, gamma, rbf, f)
h = 0.001;


xph = [Xte(:,1)+h,Xte(:,2)];
xmh = [Xte(:,1)-h,Xte(:,2)];
yph = [Xte(:,1),Xte(:,2)+h];
ymh = [Xte(:,1),Xte(:,2)-h];

Xte_eval = evaluation_matrix(rbf,gamma, Xin, Xbd, Xte);
xph_eval = evaluation_matrix(rbf,gamma, Xin, Xbd, xph);
xmh_eval = evaluation_matrix(rbf,gamma, Xin, Xbd, xmh);
yph_eval = evaluation_matrix(rbf,gamma, Xin, Xbd, yph);
ymh_eval = evaluation_matrix(rbf,gamma, Xin, Xbd, ymh);

Xte_sol = Xte_eval * alpha;
xph_sol = xph_eval * alpha;
xmh_sol = xmh_eval * alpha;
yph_sol = yph_eval * alpha;
ymh_sol = ymh_eval * alpha;

disc = (xph_sol + xmh_sol + yph_sol + ymh_sol - 4*Xte_sol)/h^2;

real = f(Xte(:,1), Xte(:,2));

error = max(abs(disc-real));
end