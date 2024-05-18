clear
clc

fun = @(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
gfun = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1); -200*(x(1)^2-x(2))];
hessian = @(x) [1200*x(1)^2-400*x(2)+2, -400*x(1); -400*x(1), 200];

x0 = [2; 1];
[x, val, k] = trust_region_method(fun, gfun, hessian, x0);

fprintf('信赖区域求解: x = [%f, %f]\n', x(1), x(2));
fprintf('最小值: f(x) = %f\n', val);
fprintf('迭代次数: %d\n', k);