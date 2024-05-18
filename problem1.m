clear
clc

fun = @(x) x(1)^2 + 25*x(2)^2;
gfun = @(x) [2*x(1); 50*x(2)];

x0=[1; 2];
[x,val,k] = gradient_descent_method(fun, gfun,x0);

fprintf('最速梯度下降法求解: x = [%f, %f]\n', x(1), x(2));
fprintf('最小值: f(x) = %f\n', val);
fprintf('迭代次数: %d\n', k);