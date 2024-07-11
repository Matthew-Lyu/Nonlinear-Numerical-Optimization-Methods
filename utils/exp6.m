clear
clc

% 数据
t = [3, 6, 9, 12, 15, 18, 21, 24];  % 时间序列
y = [57.6, 41.9, 31.0, 22.7, 16.6, 12.2, 8.9, 6.5];  % 反应物的量

x0 = [60; 0.1]; % 随机选择一个合适的起始点

% 调用 L-M 方法
[x, val, k] = L_M_least_squares_method(@Fk, @JFk, x0);

% 输出结果
fprintf('拟合参数: a = %.4f, b = %.4f\n', x(1), x(2));
fprintf('目标函数优化值: %.4f, 迭代次数: %d\n', val, k);
fprintf('拟合函数形式: y = %.2fe^{-%.2ft}\n', x(1), x(2));

% 画出数据点和拟合曲线
figure;
scatter(t, y, 'filled');
hold on;
t_fit = linspace(min(t), max(t), 100);
y_fit = x(1) * exp(-x(2) * t_fit);
plot(t_fit, y_fit, 'r-');
xlabel('时间 (t)');
ylabel('反应物的量 (y)');
title('反应速率数据拟合');
legend('数据点', '指数衰减拟合');
hold off;

function f = Fk(x)
    t0 = [3, 6, 9, 12, 15, 18, 21, 24];
    y0 = [57.6, 41.9, 31.0, 22.7, 16.6, 12.2, 8.9, 6.5];
    a = x(1);
    b = x(2);
    f = y0 - a * exp(-b * t0);
    f = f(:);
end

function jf = JFk(x)
    t0 = [3, 6, 9, 12, 15, 18, 21, 24];
    a = x(1);
    b = x(2);
    jf = [-exp(-b * t0); a * t0 .* exp(-b * t0)]';
end

