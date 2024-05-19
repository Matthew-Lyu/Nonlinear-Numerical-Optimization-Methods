%%% 利用光滑牛顿法求解信赖域子问题 %%%
function [d, val, lam, k] = trust_region_quadratic_model(gk, Bk, dta)
    % 光滑牛顿法
    % 输入：
    %   gk: 当前点的梯度
    %   Bk: 当前点的 Hessian 矩阵
    %   dta: 信赖域半径
    % 输出：
    %   d: 子问题的解
    %   val: 目标函数在 d 处的值
    %   lam: 最优拉格朗日乘子
    %   k: 迭代次数

    n = length(gk); % 计算梯度维度
    gamma = 0.05; % 初始参数
    epsilon = 1.0e-6; % 停止条件
    rho = 0.6; % 步长控制参数
    sigma = 0.2; % 用于控制步长的参数
    mu0 = 0.05; % 初始步长
    lam0 = 0.05; % 初始拉格朗日乘子
    d0 = ones(n,1); % 初始解向量
    u0 = [mu0, zeros(1,n+1)]'; % 初始辅助向量
    z0 = [mu0, lam0, d0']'; % 初始组合向量
    k = 0; % 迭代次数
    z = z0; 
    mu = mu0; 
    lam = lam0; 
    d = d0;
    
    while (k <= 150)
        % 计算子问题解的残差
        dh = dah(mu, lam, d, gk, Bk, dta);
        % 判断残差是否满足停止条件
        if (norm(dh) < epsilon)
            break;
        end
        % 构建雅可比矩阵和右侧向量
        A = JacobiH(mu, lam, d, Bk, dta);
        b = beta(mu, lam, d, gk, Bk, dta, gamma) * u0 - dh;
        % 求解线性方程组
        B = inv(A);   
        dz = B * b;
        dmu = dz(1); 
        dlam = dz(2); 
        dd = dz(3:n+2);
        m = 0;  
        mk = 0;
        % 用Armijo条件搜索合适的步长
        while (m < 20)
            dhnew = dah(mu + rho^m * dmu, lam + rho^m * dlam, d + rho^m * dd, gk, Bk, dta);
            if (norm(dhnew) <= (1 - sigma * (1 - gamma * mu0) * rho^m) * dh)
                mk = m;
                break;
            end
            m = m + 1;
        end
        alpha = rho^mk; % 更新步长
        mu = mu + alpha * dmu; % 更新步长参数
        lam = lam + alpha * dlam; % 更新拉格朗日乘子
        d = d + alpha * dd; % 更新解向量
        k = k + 1; % 迭代次数加一
    end
    val = gk' * d + 0.5 * d' * Bk * d; % 计算目标函数值
    
    function p = phi(mu, a, b)
        p = a + b - sqrt((a - b)^2 + 4 * mu); % 定义 phi 函数
    end

    function dh = dah(mu, lam, d, gk, Bk, dta)
        n1 = length(d);
        dh(1) = mu;  
        dh(2) = phi(mu, lam, dta^2 - norm(d)^2); % 计算残差中的部分
        mh = (Bk + lam * eye(n1)) * d + gk; % 计算残差中的部分
        for (i = 1:n1)
            dh(2 + i) = mh(i); % 计算残差中的部分
        end
        dh = dh(:); % 将残差向量转换为列向量
    end

    function bet = beta(mu, lam, d, gk, Bk, dta, gamma)
        dhh = dah(mu, lam, d, gk, Bk, dta); % 计算残差
        bet = gamma * norm(dhh) * min(1, norm(dhh)); % 计算步长控制参数
    end

    % 定义雅可比矩阵
    function A = JacobiH(mu, lam, d, Bk, dta)
        n2 = length(d);
        A = zeros(n2 + 2, n2 + 2); % 初始化雅可比矩阵
        pmu = -4 * mu / sqrt((lam + norm(d)^2 - dta^2)^2 + 4 * mu^2); % 计算雅可比矩阵中的元素
        thetak = (lam + norm(d)^2 - dta^2) / sqrt((lam + norm(d)^2 - dta^2)^2 + 4 * mu^2); % 计算雅可比矩阵中的元素
        A = [1,             0,            zeros(1, n2);
             pmu,          1 - thetak,  -2 * (1 + thetak) * d';
             zeros(n2, 1), d,            Bk + lam * eye(n2)]; % 构建雅可比矩阵
    end
end
