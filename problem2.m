X0 = [0; 1];
[final_x, min_val, iterations] = damp_newton_method(@fun, @gfun, @hess, X0);
disp(final_x)
disp(min_val)
disp(iterations)

X0 = [0; 1];
[final_x, min_val, iterations] = revise_newton_method(@fun, @gfun, @hess, X0);
disp(final_x)
disp(min_val)
disp(iterations)

function val = fun(X)
    x1 = X(1);
    x2 = X(2);
    val = x1^2 + x2^2 - x1*x2 - 10*x1 - 4*x2 + 60;
end

function grad = gfun(X)
    x1 = X(1);
    x2 = X(2);
    grad = [2*x1 - x2 - 10;
            2*x2 - x1 - 4];
end

function hess = hess(X)
    x1 = X(1);
    x2 = X(2);
    hess = [2, -1;
            -1, 2];
end
