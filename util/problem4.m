clear
clc

fun = @(x) x(1)^2 + x(2)^2 -1;
gfun = @(x) [2*x(1); 2*x(2)];

x = [2;2];
d = [-1;-1];
[mk, alpha, fk, newfk] = armijo_search_method(x, d, fun, gfun);
fprintf('使用Armijo方法搜索步长: alpha = %f\n', alpha);