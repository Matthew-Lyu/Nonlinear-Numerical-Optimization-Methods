clear
clc

% 参数设置
H = 2 * eye(2);
A = [1 1];
b = 1;
c = zeros(2, 1);

[x, lam, fval] = quadratic_Lagrangian(H, A, b, c);

% 检查 x2 是否满足不等式约束
if x(2) > 1/4
    % 若 x2 > 1/4，需要调整为边界情况 x2 = 1/4 并解 x1
    x(2) = 1/4;
    x(1) = 1 - x(2);
    fval = 0.5 * x' * H * x + c' * x;  % 重新计算最优值
end

% 输出结果
fprintf('最优解: (%f, %f)\n', x(1), x(2));
fprintf('最优值: %f\n', fval);