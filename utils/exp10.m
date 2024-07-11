clear
clc

% 定义目标函数
objectiveFunction = @(x) x .* sin(10 * x) + 1;

% 遗传算法参数
populationSize = 100; % 种群大小
numGenerations = 1000; % 最大迭代代数
crossoverRate = 0.8; % 交叉概率
mutationRate = 0.05; % 变异概率
tolerance = 1e-6; % 收敛精度
lowerBound = -1; % 变量下界
upperBound = 2.0; % 变量上界

% 调用遗传算法
[bestSolution, maxFitness] = genetic_algorithm(objectiveFunction, populationSize, numGenerations, crossoverRate, mutationRate, tolerance, lowerBound, upperBound);

fprintf('最优解: x = %.6f, f(x) = %.6f\n', bestSolution, maxFitness);
