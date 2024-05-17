function [x, val, k] = damp_newton_method(fun, gfun, Hess, x0)
    maxk = 100;
    rho = 0.55; sigma = 0.4;
    k = 0; epsilon = 0.01;
    while k < maxk
        gk = feval(gfun, x0); % 计算梯度
        Hk = feval(Hess, x0); % 计算海森矩阵
        dk = -inv(Hk) * gk; % 计算搜索方向
        
        m = 0; mk = 0;
        while m < 100 % 添加迭代次数限制，防止无限循环
            if feval(fun, x0 + rho^m * dk) < feval(fun, x0) + sigma * rho^m * gk' * dk
                mk = m; % 更新 mk
                break;
            end
            m = m + 1;
        end
        
        x0 = x0 + rho^mk * dk; % 更新 x0
        k = k + 1; % 更新迭代次数
        
        if norm(gk) < epsilon % 判断停止条件：梯度小于阈值
            break;
        end
    end
    x = x0;
    val = feval(fun, x);
end
