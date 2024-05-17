function [x_opt, iter] = conjugate_gradient_method(A, b, x0, epsilon, max_iter)
    x = x0;
    r = b - A * x;
    p = r;
    iter = 0;
    while true
        alpha = (r' * r) / (p' * A * p);
        x = x + alpha * p;
        r_old = r;
        r = r - alpha * A * p;
        beta = (r' * r) / (r_old' * r_old);
        p = r + beta * p;
        iter = iter + 1;
        if norm(r) < epsilon || iter >= max_iter
            break;
        end
    end
    x_opt = x;
end