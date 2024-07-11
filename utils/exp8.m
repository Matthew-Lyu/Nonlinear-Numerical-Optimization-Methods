clear
clc

% 目标函数
fun = @(x) x(1)^2 + x(2)^2;

% 梯度函数
gfun = @(x) [2 * x(1); 2 * x(2)];

% 海森矩阵函数
hessian = @(x) [2, 0; 0, 2];

% 约束函数 x1 + x2 = 2
cons = @(x) x(1) + x(2) - 2;

x0 = [10; 10]; % 初始点
tol = 1e-8; % 精度

[x_opt, f_opt, iter] = exterior_penalty_method(fun, gfun, hessian, cons, x0, tol);

fprintf('最优解: (%f, %f)\n', x_opt(1), x_opt(2));
fprintf('最优值: %f\n', f_opt);
fprintf('迭代次数: %d\n', iter);
