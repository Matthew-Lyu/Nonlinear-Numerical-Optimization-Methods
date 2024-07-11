function [x_opt, val_opt, itr] = L_M_least_squares_method(Fk, JFk, x0)
    % Levenberg-Marquardt 最小二乘法
    % 输入:
    %   x0 - 迭代的初始点
    %   Fk - 计算 F(x) 的函数句柄
    %   JFk - 计算雅可比矩阵 F'(x) 的函数句柄
    % 输出:
    %   x - F(x) 最小化时的近似解
    %   val - 解处 0.5*norm(F(x))^2 的值
    %   k - 执行的迭代次数

    % 最大迭代次数，防止无限循环
    maxk = 100;
    % Armijo 规则的参数
    rho = 0.55;  % 步长减小因子
    sigma = 0.4;  % 足够减小常数
    % 初始 mu 参数，基于残差的范数缩放
    muk = norm(feval(Fk, x0));
    % 迭代计数器
    itr = 0;
    % 梯度范数的收敛容忍度
    epsilon = 1e-6;
    % 问题的维度（变量数）
    n = length(x0);

    % 主迭代循环
    while itr < maxk
        % 在当前点计算函数值和雅可比矩阵
        fk = feval(Fk, x0);  % 当前函数值 F(x)
        jfk = feval(JFk, x0);  % 当前雅可比矩阵 F'(x)
        % 计算目标函数 0.5*norm(F(x))^2 的梯度
        gk = jfk' * fk;  % 梯度计算

        % 解算搜索方向
        dk = -(jfk' * jfk + muk * eye(n)) \ gk;  % 搜索方向

        % 检查收敛性（梯度范数是否足够小）
        if norm(gk) < epsilon, break; end

        % Armijo 线搜索确定步长
        m = 0;  % 步长指数计数器
        mk = 0;  % 找到的最佳步长指数
        while m < 20
            % 在新的候选点评估函数
            newf = 0.5 * norm(feval(Fk, x0 + rho^m * dk))^2;
            oldf = 0.5 * norm(feval(Fk, x0))^2;
            % 检查 Armijo 条件
            if newf < oldf + sigma * rho^m * gk' * dk
                mk = m;  % 接受步长
                break;
            end
            m = m + 1;
        end

        % 更新当前点
        x0 = x0 + rho^mk * dk;
        % 基于新的函数值更新 mu
        muk = norm(feval(Fk, x0));
        % 增加迭代计数器
        itr = itr + 1;
    end

    % 输出结果
    x_opt = x0;  % 最终近似解
    val_opt = 0.5 * muk^2;  % 解处的目标函数值
end
