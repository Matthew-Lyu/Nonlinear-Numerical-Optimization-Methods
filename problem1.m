x0=[1 2]';
[x,val,k] = gradient_descent(@fun,@gfun,x0); % 注意这里传递函数句柄
disp(x);
disp(val);
disp(k);

function f = fun(x)
f = x(1)^2 + 25*x(2)^2;
end

function g = gfun(x)
g = [2*x(1), 50*x(2)]';
end