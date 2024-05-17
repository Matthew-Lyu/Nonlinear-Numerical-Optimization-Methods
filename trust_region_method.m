function [x_opt, f_opt, iter] = trust_region_method(f, grad_f, hess_f, x0, delta_max, epsilon, max_iter)
    x = x0;
    delta = delta_max;
    iter = 0;
    while true
        grad = grad_f(x);
        hess = hess_f(x);
        % 解信赖域子问题
        [dx, ~] = trust_region_subproblem(grad, hess, delta);
        % 计算接受比例 rho
        rho = (f(x) - f(x + dx)) / (-grad' * dx - 0.5 * dx' * hess * dx);
        % 根据接受比例更新参数
        if rho < 0.25
            delta = 0.25 * delta;
        elseif rho > 0.75 && norm(dx) == delta
            delta = min(2 * delta, delta_max);
        end
        if rho > epsilon || norm(grad) < epsilon || iter >= max_iter
            break;
        end
        x = x + dx;
        iter = iter + 1;
    end
    x_opt = x;
    f_opt = f(x);
end