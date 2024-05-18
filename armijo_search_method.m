function [mk, alpha, fk, newfk] = armijo_search(xk, dk, fun, gfun)
    rho = 0.1;
    sigma = 0.2;
    m = 0;
    mmax = 20;
    while (m <= mmax)
        if(fun(xk + rho^m*dk) <= fun(xk) + sigma*rho^m*gfun(xk)'*dk)
            mk = m; break;
        end
        m = m+1;
    end
    alpha = rho^mk;
    newxk = xk + alpha*dk;
    fk = fun(xk);
    newfk = fun(newxk);
end
