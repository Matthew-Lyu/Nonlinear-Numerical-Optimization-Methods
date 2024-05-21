function [mk, alpha, fk, newfk] = armijo_search_method(xk, dk, fun, grad)
    % Armijo准则搜索
    % 输入：
    %   xk: 当前点
    %   dk: 搜索方向
    %   fun: 目标函数
    %   grad: 梯度函数
    % 输出：
    %   mk: 步长指数
    %   alpha: 步长因子
    %   fk: 当前点的函数值
    %   newfk: 新点的函数值

    rho = 0.1; % 步长缩放因子
    sigma = 0.2; % Armijo条件系数
    m = 0; % 步长指数
    mmax = 20; % 最大步长指数

    while (m <= mmax)
        if fun(xk + rho^m * dk) <= fun(xk) + sigma * rho^m * grad(xk)' * dk
            mk = m; % 更新步长指数
            break;
        end
        m = m + 1;
    end

    alpha = rho^mk; % 计算步长因子
    newxk = xk + alpha * dk; % 计算新点
    fk = fun(xk); % 当前点的函数值
    newfk = fun(newxk); % 新点的函数值
end
