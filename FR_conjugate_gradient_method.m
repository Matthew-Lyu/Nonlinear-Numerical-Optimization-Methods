function [x_opt, f_opt, iter] = FR_conjugate_gradient_method(fun, grad, x0)
    % FR共轭梯度法
    % 输入：
    %   fun: 目标函数
    %   grad: 梯度函数
    %   x0: 初始点
    % 输出：
    %   x_opt: 近似最优解
    %   f_opt: 最优解对应的函数值
    %   iter: 迭代次数

    max_iter = 1000; % 最大迭代次数
    rho = 0.6; % 线搜索参数
    sigma = 0.4; % Armijo条件系数
    epsilon = 0.01; % 误差阈值

    n = length(x0);
    k = 0; % 迭代次数

    while k < max_iter
        % 计算当前点的梯度
        g = feval(grad, x0);

        % 计算FR系数
        itern = k - (n + 1) * floor(k / (n + 1)) + 1;

        % 计算搜索方向
        if itern == 1
            d = -g;
        else
            beta = (g' * g) / (g0' * g0);
            d = -g + beta * d0;
            gd = g' * d;
            if gd >= 0.0
                d = -g;
            end
        end

        % 检验终止条件：梯度小于阈值
        if norm(g) < epsilon
            break;
        end

        % Armijo搜索
        m = 0;
        mk = 0;
        while m < 20
            if feval(fun, x0 + rho^m * d) < feval(fun, x0) + sigma * rho^m * g' * d
                mk = m;
                break;
            end
            m = m + 1;
        end

        % 更新当前点和迭代次数
        x0 = x0 + rho^mk * d;
        f_opt = feval(fun, x0);
        g0 = g;
        d0 = d;
        k = k + 1;
    end

    % 输出最优解和最优值
    x_opt = x0;
    f_opt = feval(fun, x0);
    iter = k;
end
