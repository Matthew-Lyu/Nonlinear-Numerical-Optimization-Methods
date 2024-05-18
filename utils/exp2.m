clear
clc

fun = @(x) x(1)^2 + x(2)^2 - x(1)*x(2)- 10*x(1) - 4*x(2) + 60;
gfun = @(x) [2*x(1) - x(2) - 10; 2*x(2) - x(1) - 4];
hessian = @(X) [2, -1; -1, 2];

x0 = [1; 2];
[x,val,k] = damp_newton_method(fun, gfun, hessian, x0);
fprintf('阻尼牛顿法求解: x = [%f, %f]\n', x(1), x(2));
fprintf('最小值: f(x) = %f\n', val);
fprintf('迭代次数: %d\n', k);

disp(' ')

[x,val,k] = revised_newton_method(fun, gfun, hessian, x0);
fprintf('修正牛顿法求解: x = [%f, %f]\n', x(1), x(2));
fprintf('最小值: f(x) = %f\n', val);
fprintf('迭代次数: %d\n', k);


