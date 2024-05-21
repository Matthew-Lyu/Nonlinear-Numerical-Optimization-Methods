function [x_opt, f_opt, iter] = gradient_descent_method(fun, grad_fun, x0)
    % 最速梯度下降法
    % 输入：
    %   fun: 目标函数
    %   grad_fun: 梯度函数
    %   x0: 初始点
    % 输出：
    %   x_opt: 近似最优解
    %   f_opt: 最优解对应的函数值
    %   iter: 迭代次数

    max_iter = 5000;  % 最大迭代次数
    rho = 0.5;  % Armijo线搜索参数
    sigma = 0.4;  % Armijo线搜索参数
    iter = 0;  % 迭代次数
    epsilon = 0.01;  % 误差阈值

    % 初始化
    x = x0;

    while(iter < max_iter)
        % 计算梯度
        g = feval(grad_fun, x);

        % 计算搜索方向为负梯度方向
        d = -g;

        % 若梯度范数小于误差阈值，则停止迭代
        if(norm(d) < epsilon)
            break;
        end

        % 初始化线搜索参数
        m = 0;
        mk = 0;

        % 线搜索确定步长因子
        while(m < 20)
            % 判断是否满足Armijo条件
            if(feval(fun, x + rho^m * d) < feval(fun, x) + sigma * rho^m * g' * d)
                mk = m;
                break;
            end
            m = m + 1;
        end

        % 更新当前点
        x = x + rho^mk * d;
        iter = iter + 1; 
    end
    
    % 输出最优解和最优值
    x_opt = x;
    f_opt = feval(fun, x);
end

