function [minError, errorTrace] = CoDE_Test(NP, D, maxFES, searchRange, fhd, funcNum, realMinVal)
%CODE_TEST CoDE算法测试函数，用于测试算法性能。
% 输入：
% NP：种群数量，D：维度，maxFES：最大函数评估次数，searchRange：搜索范围（1*2），
% fhd：测试函数句柄，funcNum：测试函数序号，realMinVal：真实最小值。
% 输出：
% [minError, errorTrace]
% minError：最小误差值，
% errorTrace：每一代最小误差值记录（1*(G+1))。


% 参数库
F = [1, 1, 0.8];
CR = [0.1, 0.9, 0.2];

% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
xCost = fhd(x, funcNum);  % 初始成本
FES = NP;
u = zeros(D, 3);  % 3列分别储存3种策略生成的试验个体

G = ceil((maxFES - NP) / NP / 3);  % 最大代数预分配内存
errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = abs(min(xCost) - realMinVal);

% 迭代
g = 1;
while FES < maxFES
    for i = 1 : NP
        % DE/rand/1/bin
        parameterIndex = randi(3);  % 选择参数组合
        r1 = randi(NP);
        while r1 == i
            r1 = randi(NP);
        end
        r2 = randi(NP);
        while r2 == i || r2 == r1
            r2 = randi(NP);
        end
        r3 = randi(NP);
        while r3 == i || r3 == r1 || r3 == r2
            r3 = randi(NP);
        end
        jRand = randi(D);
        % 变异交叉
        for j = 1 : D
            if rand() <= CR(parameterIndex) || j == jRand
                u(j, 1) = x(j, r1) + F(parameterIndex) * (x(j, r2) - x(j, r3));
            else
                u(j, 1) = x(j, i);
            end
        end

        % DE/rand/2/bin
        parameterIndex = randi(3);  % 选择参数组合
        r1 = randi(NP);
        while r1 == i
            r1 = randi(NP);
        end
        r2 = randi(NP);
        while r2 == i || r2 == r1
            r2 = randi(NP);
        end
        r3 = randi(NP);
        while r3 == i || r3 == r1 || r3 == r2
            r3 = randi(NP);
        end
        r4 = randi(NP);
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3
            r4 = randi(NP);
        end
        r5 = randi(NP);
        while r5 == i || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4
            r5 = randi(NP);
        end
        jRand = randi(D);
        % 变异交叉
        for j = 1 : D
            if rand() <= CR(parameterIndex) || j == jRand
                u(j, 2) = x(j, r1) + rand() * (x(j, r2) - x(j, r3)) + F(parameterIndex) * (x(j, r4) - x(j, r5));  % 第一个F使用[0,1]随机数提高搜索能力
            else
                u(j, 2) = x(j, i);
            end
        end

        % DE/current-to-rand/1
        parameterIndex = randi(3);  % 选择参数组合
        r1 = randi(NP);
        while r1 == i
            r1 = randi(NP);
        end
        r2 = randi(NP);
        while r2 == i || r2 == r1
            r2 = randi(NP);
        end
        r3 = randi(NP);
        while r3 == i || r3 == r1 || r3 == r2
            r3 = randi(NP);
        end
        u(:, 3) = x(:, i) + rand() .* (x(:, r1) - x(:, i)) + F(parameterIndex) .* (x(:, r2) - x(:, r3));
        
        % 越界截断
        u(u < searchRange(1)) = searchRange(1);
        u(u > searchRange(2)) = searchRange(2);

        % 选择
        if maxFES - FES < 3
            % 剩余评估次数不足3次时
            uCost = fhd(u(:, 1 : maxFES - FES), funcNum);
            FES = maxFES;
        else
            uCost = fhd(u, funcNum);
            FES = FES + 3;
        end
        [minUCost, index] = min(uCost);  % 取最优策略实验个体
        if minUCost <= xCost(i)
            x(:, i) = u(:, index);
            xCost(i) = minUCost;
        end

        % 中途强行跳出
        if FES == maxFES
            break
        end
    end
    errorTrace(g + 1) = abs(min(xCost) - realMinVal);
    g = g + 1;
end

minError = errorTrace(end);

end
