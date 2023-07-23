function [minVal, solution, trace] = CoDE(NP, D, G, searchRange, fhd, funcNum)
%CODE CoDE算法函数。
% 输入：
% NP：种群数量，D：维度，G：代数，searchRange：搜索范围（1*2），
% fhd：测试函数句柄，funcNum：测试函数序号。
% 输出：
% [minVal, solution, trace]
% minVal：函数最小值，
% solution：函数最优解（D*1）,
% trace：每一代的函数最小值（1*(G+1))。


% 参数库
F = [1, 1, 0.8];
CR = [0.1, 0.9, 0.2];

% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
xCost = fhd(x, funcNum);  % 初始成本
u = zeros(D, 3);  % 3列分别储存3种策略生成的试验个体


trace = zeros(1, G + 1);  % 储存每代最小值
trace(1) = min(xCost);

% 迭代
for g = 1 : G
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
                % 越界调整
                if u(j, 1) < searchRange(1)
                    u(j, 1) = (searchRange(1) + x(j, i)) / 2;
                elseif u(j, 1) > searchRange(2)
                    u(j, 1) = (searchRange(2) + x(j, i)) / 2;
                end
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
                % 越界调整
                if u(j, 2) < searchRange(1)
                    u(j, 2) = (searchRange(1) + x(j, i)) / 2;
                elseif u(j, 2) > searchRange(2)
                    u(j, 2) = (searchRange(2) + x(j, i)) / 2;
                end
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
        % 越界调整
        for j = 1 : D
            if u(j, 3) < searchRange(1)
                u(j, 3) = (searchRange(1) + x(j, i)) / 2;
            elseif u(j, 3) > searchRange(2)
                u(j, 3) = (searchRange(2) + x(j, i)) / 2;
            end
        end

        % 选择
        uCost = fhd(u, funcNum);
        [minUCost, index] = min(uCost);  % 取最优策略实验个体
        if minUCost <= xCost(i)
            x(:, i) = u(:, index);
            xCost(i) = minUCost;
        end
    end
    trace(g + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);

end

