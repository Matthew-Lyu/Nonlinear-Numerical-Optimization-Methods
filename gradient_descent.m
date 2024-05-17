function [x, num_iterations] = gradient_descent(x0, epsilon)
    % 初始化参数
    x = x0;
    num_iterations = 0;
    
    % 计算初始梯度
    grad = [2*x(1); 50*x(2)];
    
    % 迭代直到梯度的模长小于 epsilon
    while norm(grad) > epsilon
        % 计算步长，这里使用一个简单的步长选择
        alpha = 1 / (1 + num_iterations);
        
        % 更新位置
        x = x - alpha * grad;
        
        % 更新梯度
        grad = [2*x(1); 50*x(2)];
        
        % 迭代次数加一
        num_iterations = num_iterations + 1;
    end
end
