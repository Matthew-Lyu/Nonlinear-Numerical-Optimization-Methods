function gradient_descent_3d_gif()
    % 定义目标函数和梯度
    f = @(x, y) x.^2 + y.^2; % 函数
    grad_f = @(x, y) [2*x; 2*y]; % 梯度

    % 初始化
    x0 = [3; 2]; % 初始点
    alpha = 0.1; % 步长
    num_iterations = 20; % 迭代次数

    % 创建图形
    figure;
    hold on;
    axis tight manual;
    zlim([0, 3]);
    view(-60, 30);

    % 生成网格数据用于绘图
    [X, Y] = meshgrid(-4:0.1:4, -3:0.1:3);
    Z = f(X, Y);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    colormap jet;

    % GIF文件设置
    filename = 'gradient_descent_v2.gif';

    % 迭代优化过程
    x = x0;
    for i = 1:num_iterations
        % 计算函数值和梯度
        fx = f(x(1), x(2));
        grad = grad_f(x(1), x(2));

        % 绘制当前点
        plot3(x(1), x(2), fx, 'ko', 'MarkerFaceColor', 'k');
        % 绘制梯度方向
        quiver3(x(1), x(2), fx, -alpha*grad(1), -alpha*grad(2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 3);
        title(['Iteration ', num2str(i)]);
        drawnow;

        % 更新位置
        x = x - alpha * grad;

        % 捕获图形并写入GIF
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
        end
    end

    hold off;
end
