function [x_optimal, optimal_value, iterations] = trust_region_method(fun, gfun, Hess, x_initial)
    % 信赖域方法求解无约束优化问题
    % 输入：
    %   fun: 目标函数
    %   gfun: 梯度函数
    %   Hess: Hessian矩阵函数
    %   x_initial: 初始点
    % 输出：
    %   x_optimal: 近似最优解
    %   optimal_value: 目标函数在最优解处的值
    %   iterations: 迭代次数

    % 设置参数
    eta1 = 0.1;
    eta2 = 0.75;
    dtabar = 2.0;
    tau1 = 0.5;
    tau2 = 2.0;
    epsilon = 1e-6;

    % 初始化
    x = x_initial;
    delta = 1;
    Bk = feval(Hess, x);

    iterations = 0;
    while iterations < 100000
        gk = feval(gfun, x);
        
        % 检查梯度是否满足终止条件
        if norm(gk) < epsilon
            break;
        end

        % 解信赖域子问题
        [d, val, lambda, ik] = trust_region_quadratic_model(gk, Bk, delta);

        % 计算实际下降量和预测下降量
        delta_q = -qk(x, d);
        delta_f = fun(x) - fun(x + d);
        rk = delta_f / delta_q;

        % 更新信赖域半径
        if rk <= eta1
            delta = tau1 * delta;
        else 
            if rk >= eta2 && norm(d) == delta
                delta = min(tau2 * delta, dtabar);
            end
        end

        % 更新迭代点
        if rk > eta1
            x = x + d;
            Bk = feval(Hess, x);
        end

        iterations = iterations + 1; % 更新迭代次数
    end

    x_optimal = x;
    optimal_value = fun(x_optimal);

    % 信赖域子问题目标函数
    function q_d = qk(x, d)
        g_k = feval(gfun, x);
        B_k = feval(Hess, x);
        q_d = g_k' * d + 0.5 * d' * B_k * d;
    end
end
