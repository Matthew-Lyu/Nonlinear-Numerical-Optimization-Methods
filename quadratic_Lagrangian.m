function [x_opt, lam, f_opt] = quadratic_Lagrangian(H, A, b, c)
    % 拉格朗日二次规划函数
    % 用于求解具有等式约束的二次优化问题。
    % 形式为: min f(x) = 0.5*x'*H*x + c'*x, 其中约束条件为 Ax = b。
    %
    % 输入:
    % H - 目标函数中二次项的系数矩阵，应为正定矩阵。
    % A - 约束条件中的系数矩阵。
    % b - 约束条件中的常数向量。
    % c - 目标函数中一次项的系数向量。
    %
    % 输出:
    % x - 最优解向量。
    % lam - 拉格朗日乘子向量。
    % f_opt - 在x处的目标函数最小值。
    
    % 计算 H 的逆矩阵
    inv_H = inv(H);
    
    % 计算 A*H^-1*A'，即约束矩阵 A，预乘 H 的逆，再后乘 A 的转置
    AHA = A * inv_H * A';
    
    % 计算 A*H^-1*A' 的逆矩阵
    inv_AHA = inv(AHA);
    
    % 计算 A*H^-1
    AiH = A * inv_H;
    
    % 计算 G 矩阵，用于后续计算最优解向量 x
    G = inv_H - AiH' * inv_AHA * AiH;
    
    % 计算 B 矩阵，用于后续计算拉格朗日乘子 lam
    B = inv_AHA * AiH;
    
    % 计算 C 矩阵，也是用于后续计算拉格朗日乘子 lam
    C = -inv_AHA;
    
    % 计算最优解向量 x
    x_opt = B' * b - G * c;
    
    % 计算拉格朗日乘子向量 lam
    lam = B * c - C * b;
    
    % 计算在 x 处的目标函数值
    f_opt = 0.5 * x_opt' * H * x_opt + c' * x_opt;

end
