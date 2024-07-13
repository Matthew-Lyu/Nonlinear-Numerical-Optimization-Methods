function [x_opt, f_opt, iter] = interior_penalty_method(fun, grad, hess, cons,cons_grad,cons_hessian, x0, tol)
    % 内点罚函数法
    % 该方法通过在目标函数中加入障碍项（惩罚项）来确保解不越过约束界限。
    % 输入：
    %   fun - 原始的目标函数
    %   grad - 原始目标函数的梯度
    %   hess - 原始目标函数的海森矩阵
    %   cons - 约束函数
    %   cons_grad - 约束函数的梯度
    %   cons_hessian - 约束函数的海森矩阵
    %   x0 - 优化问题的初始点
    %   tol - 容忍度，用于控制优化的精度
    % 输出：
    %   x_opt - 最优解
    %   f_opt - 在最优解处的目标函数值
    %   iter - 实际迭代次数

    R = 10;  % 初始障碍因子，用于控制障碍项的强度
    c = 0.7; % 障碍因子的递减系数，用于每轮迭代后减少障碍因子
    x = x0;  % 设置优化算法的起始点为输入的初始点
    max_iter = 100; % 设置最大迭代次数
    iter = 0;

    while iter < max_iter
        % 定义内点罚函数的障碍项，适用于输入的约束函数
        barrier = @(x) -R * log(cons(x));
        % 定义内点罚函数的梯度
        barrier_grad = @(x) -R ./ cons(x) .* cons_grad(x);
        % 定义内点罚函数的海森矩阵
        barrier_hessian = @(x) R * (cons_grad(x) * cons_grad(x)' ./ cons(x).^2 - cons_hessian(x) ./ cons(x));

        % 定义内点罚函数，结合原始目标函数和障碍项
        penalty_fun = @(x) fun(x) + barrier(x);
        % 定义内点罚函数的梯度
        penalty_grad = @(x) grad(x) + barrier_grad(x);
        % 定义内点罚函数的海森矩阵
        penalty_hess = @(x) hess(x) + barrier_hessian(x);

        % 使用阻尼牛顿法求解当前无约束优化问题
        [x, ~, ~] = damp_newton_method(penalty_fun, penalty_grad, penalty_hess, x);

        % 检查障碍项的绝对值是否已小于容忍度，以判断是否接近约束边界
        if abs(barrier(x)) < tol
            break; % 若满足容忍度，则终止迭代
        else
            R = R * c; % 若不满足，递减障碍因子，继续迭代
        end

        iter = iter + 1; % 更新迭代次数
    end

    % 计算在当前解x处的目标函数值
    f_opt = fun(x);
    % 返回最优解
    x_opt = x;
end