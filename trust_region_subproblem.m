function [dx, lambda] = trust_region_subproblem(grad, hess, delta)
    n = length(grad);
    % 构造二次子问题
    A = hess;
    b = grad;
    H = [A, -eye(n); -eye(n), zeros(n)];
    f = zeros(2 * n, 1);
    % 设置线性约束
    lb = -delta * ones(n, 1);
    ub = delta * ones(n, 1);
    % 使用 MATLAB 自带的二次规划求解器 quadprog
    [dx, ~, ~, ~, lambda] = quadprog(H, f, [], [], [], [], lb, ub, []);
    dx = dx(1:n);
end