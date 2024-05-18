function [x,val,k] = gradient_descent(fun,gfun,x0)

%功能: 用最速梯度下降法求解无约束问题: min f(x)
%输入: x0是初始点, fun, gfun分别是目标函数和梯度 %输出: x, val分别是近似最优点和最优值, k是迭代次数. 

maxk=5000; %最大迭代次数
k=0; epsilon=0.01;
while(k < maxk)
    g=feval(gfun,x0); %计算梯度 
    d=-g; %计算搜索方向 
    if(norm(d) < epsilon), break; end 
    
    % 精确线搜索（一维搜索）
    alpha = exact_line_search(fun, gfun, x0, d);
    
    x0 = x0 + alpha * d;
    k=k+1; 
end
x=x0;
val=feval(fun,x0);

function alpha = exact_line_search(fun, gfun, x0, d)
    alpha = 1; % 初始步长
    c = 0.5; % Armijo搜索参数
    rho = 0.5; % 步长衰减因子
    
    while true
        if feval(fun, x0 + alpha * d) <= feval(fun, x0) + c * alpha * gfun(x0)' * d
            break;
        else
            alpha = rho * alpha; % 步长衰减
        end
    end
end

end
