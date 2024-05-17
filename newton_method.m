function [x_opt, iter] = newton_method(f, grad_f, hess_f, x0, epsilon, max_iter)
    x = x0;
    iter = 0;
    while true
        grad = grad_f(x);
        hess = hess_f(x);
        dx = -inv(hess) * grad;
        x = x + dx;
        iter = iter + 1;
        if norm(grad) < epsilon || iter >= max_iter
            break;
        end
    end
    x_opt = x;
end