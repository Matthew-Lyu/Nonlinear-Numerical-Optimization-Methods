clear
clc

% 目标函数
fun = @(x) x(1)^2 + x(2)^2;

% 目标函数的梯度
gfun = @(x) [2 * x(1); 2 * x(2)];

% 目标函数的海森矩阵
hessian = @(x) [2, 0; 0, 2];

% 内点罚函数的障碍项，适用于 x1 - 1 >= 0
barrier = @(x, R) -R * log(x(1) - 1);

% 约束的梯度
barrier_grad = @(x, R) [-R / (x(1) - 1); 0];

% 约束的海森矩阵
barrier_hessian = @(x, R) [R / (x(1) - 1)^2, 0; 0, 0];

x0 = [10; 10]; % 初始点
tol = 1e-8; % 精度要求

[x_opt, f_opt, iter] = interior_penalty_method(fun, barrier, gfun, hessian,barrier_grad,barrier_hessian, x0, tol);

fprintf('最优解: (%f, %f)\n', x_opt(1), x_opt(2));
fprintf('最优值: %f\n', f_opt);
fprintf('迭代次数: %d\n', iter);
