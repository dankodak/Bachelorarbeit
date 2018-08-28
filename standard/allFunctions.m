function [rbf, lap_rbf, lap2_rbf, f, w, realSol, realSolPlot] = allFunctions(kernel, pde, symmetric)
syms a b c d gammas rbfs lap_rbfs ws prods fs realSols realSolPlots
%% PDE
switch pde
    case 'square'
        ws(a,b)= -a^2 - b^2 + 2 -sqrt(a^4 - 2*a^2 + b^4 - 2*b^2 +2);
        fs(a,b) = - 2*pi^2*sin(pi*a)*sin(pi*b);
%         realSols(a,b) = sin(pi*a)*sin(pi*b);
        realSols(a,b) = 0;
        realSolPlots(a,b) = sin(pi*a)*sin(pi*b);
        
        w = matlabFunction(ws);
        f = matlabFunction(fs);
        realSol = matlabFunction(realSols);
        realSolPlot = matlabFunction(realSolPlots);
    case 'circle'
        ws(a,b) = 1 - a^2 - b^2;
        fs(a,b) = exp(-a^2-b^2)*(-4+4*(a^2+b^2));
        % realSols(a,b) = exp(-a^2-b^2) - 1/exp(1);
        realSols(a,b) = 0;
        realSolPlots(a,b) = exp(-a^2-b^2) - 1/exp(1);
        
        w = matlabFunction(ws);
        f = matlabFunction(fs);
        realSol = matlabFunction(realSols);
        realSolPlot = matlabFunction(realSolPlots);
end
%% Kernel
switch kernel
    case 'wendland'
        switch symmetric
            case 0  
                rbfs(gammas,a,b,c,d) = (1 - gammas*((c-a)^2 + (d-b)^2))^4 * (4*gammas*((c-a)^2 + (d-b)^2) + 1)/(20*gammas^2);
%                 rbfs(gammas,a,b,c,d) = ((1-gammas*((c-a)^2 + (d-b)^2))^6 * (3+gammas*((c-a)^2 + (d-b)^2) * (18 + 35* gammas * ((c-a)^2 + (d-b)^2))))/(1680*gammas^4);
                rbf1 = matlabFunction(rbfs);
                lap_rbf1 = matlabFunction(diff(rbfs,a,2) + diff(rbfs,b,2));
                lap2_rbf1 = 0;
                [rbf, lap_rbf, lap2_rbf] = wendland(rbf1, lap_rbf1, lap2_rbf1);
            case 1
                rbfs(gammas,a,b,c,d) =((1-gammas*((c-a)^2 + (d-b)^2))^6 * (3+gammas*((c-a)^2 + (d-b)^2) * (18 + 35* gammas * ((c-a)^2 + (d-b)^2))))/(1680*gammas^4);
                rbf1 = matlabFunction(rbfs);
                lap_rbf1s = diff(rbfs,a,2) + diff(rbfs,b,2);
                lap_rbf1 = matlabFunction(lap_rbf1s);
                lap2_rbf1 = matlabFunction(diff(lap_rbf1s,c,2) + diff(lap_rbf1s,d,2));
                [rbf, lap_rbf, lap2_rbf] = wendland(rbf1, lap_rbf1, lap2_rbf1);
                
        end
    case 'gauss'
        switch symmetric
            case 0
                rbfs(gammas,a,b,c,d) = exp(-gammas*((c-a)^2 + (d-b)^2));
                rbf = matlabFunction(rbfs);
                lap_rbf = matlabFunction(diff(rbfs,a,2) + diff(rbfs,b,2));
                lap2_rbf = 0;
            case 1
                rbfs(gammas,a,b,c,d) = exp(-gammas*((c-a)^2 + (d-b)^2));
                rbf = matlabFunction(rbfs);
                lap_rbfs = diff(rbfs,a,2) + diff(rbfs,b,2);
                lap_rbf = matlabFunction(lap_rbfs);
                lap2_rbf = matlabFunction(diff(lap_rbfs,c,2) + diff(lap_rbfs,d,2));
        end
end




%%
% rbf = matlabFunction(rbfs);
% w = matlabFunction(ws);
% f = matlabFunction(fs);
% realSol = matlabFunction(realSols);
% realSolPlot = matlabFunction(realSolPlots);
% prods(gammas,a,b,c,d) = rbfs(gammas,a,b,c,d) * w(a,b);
% lap_rbf = matlabFunction(diff(prods,a,2) + diff(prods,b,2));
% lap2_rbf = matlabFunction(diff(prods,a,4) + diff(prods,b,4) + 2*diff(diff(prods,a,2),b,2));
end