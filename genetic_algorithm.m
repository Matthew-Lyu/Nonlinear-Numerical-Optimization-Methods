function [x_opt, maxFitness] = genetic_algorithm(objectiveFunction, populationSize, numGenerations, crossoverRate, mutationRate, tolerance, lowerBound, upperBound)
    % 遗传算法
    % 输入：
    %   objectiveFunction - 目标函数，用于计算种群的适应度值。
    %   populationSize - 种群大小，即种群中个体的数量。
    %   numGenerations - 迭代次数，算法运行的代数。
    %   crossoverRate - 交叉率，决定多少比例的个体会进行交叉。
    %   mutationRate - 变异率，决定多少比例的个体会发生变异。
    %   tolerance - 收敛容忍度，当种群的适应度变化小于此值时算法停止。
    %   lowerBound - 搜索空间的下界。
    %   upperBound - 搜索空间的上界。
    % 输出：
    %   x_opt - 最优解，适应度最高的个体。
    %   maxFitness - 最大适应度值。

    % 设置二进制编码长度
    numBits = 16;
    
    % 初始化种群：在搜索空间内随机生成初始二进制种群
    population = randi([0, 1], populationSize, numBits);
    fitness = calculate_fitness(objectiveFunction, population, lowerBound, upperBound, numBits); % 计算每个个体的适应度

    % 初始化最优个体
    [maxFitness, bestIndex] = max(fitness);
    bestIndividual = population(bestIndex, :);

    for generation = 1:numGenerations
        % 选择操作（最优保存策略）：保留当前最优个体
        newPopulation = population;
        newPopulation(1, :) = bestIndividual;

        % 交叉操作：随机选择个体对进行基因交换，产生新的个体
        for i = 2:2:populationSize-1
            if rand < crossoverRate
                crossoverPoint = randi([1, numBits-1]);
                newPopulation(i, crossoverPoint+1:end) = population(i+1, crossoverPoint+1:end);
                newPopulation(i+1, crossoverPoint+1:end) = population(i, crossoverPoint+1:end);
            end
        end

        % 变异操作：随机改变个体的某些基因，以增加种群的多样性
        for i = 2:populationSize
            for j = 1:numBits
                if rand < mutationRate
                    newPopulation(i, j) = 1 - newPopulation(i, j); % 翻转基因位
                end
            end
        end

        % 更新种群和适应度：用新种群替换旧种群，并重新计算适应度
        population = newPopulation;
        fitness = calculate_fitness(objectiveFunction, population, lowerBound, upperBound, numBits);

        % 更新最优个体
        [currentMaxFitness, bestIndex] = max(fitness);
        if currentMaxFitness > maxFitness
            maxFitness = currentMaxFitness;
            bestIndividual = population(bestIndex, :);
        end

        % 打印当前代数和最大适应度
        fprintf('迭代数 %d: 当前适应度 = %.6f\n', generation, maxFitness);

        % 检查收敛：如果最大和最小适应度之差小于容忍度，则停止迭代
        if max(fitness) - min(fitness) < tolerance
            break;
        end
    end

    % 找到最优解：找出适应度最高的个体及其适应度值
    x_opt = decode_individual(bestIndividual, lowerBound, upperBound, numBits);

    fprintf('最优解: x = %.6f, f(x) = %.6f\n', x_opt, maxFitness);
end

% 计算种群的适应度值
function fitness = calculate_fitness(objectiveFunction, population, lowerBound, upperBound, numBits)
    populationSize = size(population, 1);
    fitness = zeros(populationSize, 1);
    for i = 1:populationSize
        x = decode_individual(population(i, :), lowerBound, upperBound, numBits);
        fitness(i) = objectiveFunction(x);
    end
end

% 将二进制个体解码为连续值
function x = decode_individual(individual, lowerBound, upperBound, numBits)
    binaryString = num2str(individual);
    binaryValue = bin2dec(binaryString);
    x = lowerBound + (upperBound - lowerBound) * (binaryValue / (2^numBits - 1));
end
