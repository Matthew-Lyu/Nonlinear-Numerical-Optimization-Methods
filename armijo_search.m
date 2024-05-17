function alpha = armijo_search(f, grad_f, x, d, rho, c)
    alpha = 1; % 初始步长
    while f(x + alpha * d) > f(x) + c * alpha * grad_f(x)' * d
        alpha = rho * alpha; % 减小步长
    end
end