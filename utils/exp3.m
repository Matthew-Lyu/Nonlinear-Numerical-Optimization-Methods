clear
clc

fun = @(x) x(1)^2 + 4 * x(2)^2;
gfun = @(x) [2 * x(1); 8 * x(2)];

x0 = [1; 1];
[x, val, k] = FR_conjugate_gradient_method(fun, gfun, x0);
fprintf('FR共轭梯度法求解: x = [%f, %f]\n', x(1), x(2));
fprintf('最小值: f(x) = %f\n', val);
fprintf('迭代次数: %d\n', k);