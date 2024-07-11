function [x_opt, f_val, iter] = exterior_penalty_method(fun, gfun, hfun, cons, x0, tol)
    % 外罚函数法用于解决有约束的优化问题。
    % 利用惩罚因子将约束条件整合到目标函数中，转换成无约束优化问题解决。
    % 输入：
    %   fun - 目标函数
    %   gfun - 目标函数的梯度
    %   hfun - 目标函数的海森矩阵
    %   cons - 约束函数
    %   x0 - 初始迭代点
    %   tol - 约束满足的容忍度
    % 输出：
    %   x_opt - 最优解
    %   f_val - 在x_opt处的目标函数值
    %   iter - 迭代次数

    M = 1;  % 初始惩罚因子，用于控制约束违反的惩罚程度
    c = 10; % 惩罚因子的增长率，用于调整每轮迭代后惩罚因子的大小
    x = x0; % 初始点设置为输入的初始迭代点
    max_iter = 100; % 设置最大迭代次数以防止无限循环
    iter = 0;

    % 主循环
    while iter < max_iter
        % 定义外罚函数，结合原始目标函数与约束违反惩罚
        penalty_fun = @(x) fun(x) + M * cons(x)^2;
        
        % 定义外罚函数的梯度，包括原始梯度与约束违反的梯度部分
        penalty_grad = @(x) gfun(x) + 2 * M * cons(x) * [1; 1];
        
        % 定义外罚函数的海森矩阵，增加对应的二次项以处理线性约束
        penalty_hess = @(x) hfun(x) + 2 * M * [1, 1; 1, 1];

        % 调用阻尼牛顿法来求解无约束优化问题
        [x, ~, ~] = damp_newton_method(penalty_fun, penalty_grad, penalty_hess, x);

        % 检查约束是否满足到指定的容忍度
        if abs(cons(x)) < tol
            break; % 如果满足，则退出循环
        else
            M = M * c; % 否则增加惩罚因子，更加严格地处理约束违反
        end

        iter = iter + 1; %
    end

    % 计算最终的函数值并返回解
    f_val = fun(x);
    x_opt = x;
end
