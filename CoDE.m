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
u = zeros(D, NP);
uTemp = zeros(D, 3);  % 暂存三个试验个体
xCost = fhd(x, funcNum);  % 初始成本
uCost = zeros(1, NP);


trace = zeros(1, G + 1);  % 储存每代最小值
trace(1) = min(xCost);

% 迭代
for g = 1 : G
    for i = 1 : NP
        % DE/rand/1/bin
        parameterIndex = randi(3);  % 选择参数组合
        r = randperm(NP, 3);
        regenIndex = r == i;  % 只可能有1个等于i
        while r(regenIndex) == i
            r(regenIndex) = randi(NP);
        end
        jRand = randi(D);
        % 变异交叉
        for j = 1 : D
            if rand() <= CR(parameterIndex) || j == jRand
                uTemp(j, 1) = x(j, r(1)) + F(parameterIndex) * (x(j, r(2)) - x(j, r(3)));
                % 越界调整
                if uTemp(j, 1) < searchRange(1)
                    uTemp(j, 1) = (searchRange(1) + x(j, i)) / 2;
                elseif uTemp(j, 1) > searchRange(2)
                    uTemp(j, 1) = (searchRange(2) + x(j, i)) / 2;
                end
            else
                uTemp(j, 1) = x(j, i);
            end
        end

        % DE/rand/2/bin
        parameterIndex = randi(3);  % 选择参数组合
        r = randperm(NP, 5);
        regenIndex = r == i;  % 只可能有1个等于i
        while r(regenIndex) == i
            r(regenIndex) = randi(NP);
        end
        jRand = randi(D);
        % 变异交叉
        for j = 1 : D
            if rand() <= CR(parameterIndex) || j == jRand
                uTemp(j, 2) = x(j, r(1)) + rand() * (x(j, r(2)) - x(j, r(3))) + F(parameterIndex) * (x(j, r(4)) - x(j, r(5)));  % 第一个F使用[0,1]随机数提高搜索能力
                % 越界调整
                if uTemp(j, 2) < searchRange(1)
                    uTemp(j, 2) = (searchRange(1) + x(j, i)) / 2;
                elseif uTemp(j, 2) > searchRange(2)
                    uTemp(j, 2) = (searchRange(2) + x(j, i)) / 2;
                end
            else
                uTemp(j, 2) = x(j, i);
            end
        end

        % DE/current-to-rand/1
        parameterIndex = randi(3);  % 选择参数组合
        r = randperm(NP, 3);
        regenIndex = r == i;  % 只可能有1个等于i
        while r(regenIndex) == i
            r(regenIndex) = randi(NP);
        end
        uTemp(:, 3) = x(:, i) + rand() .* (x(:, r(1)) - x(:, i)) + F(parameterIndex) .* (x(:, r(2)) - x(:, r(3)));
        % 越界调整
        for j = 1 : D
            if uTemp(j, 3) < searchRange(1)
                uTemp(j, 3) = (searchRange(1) + x(j, i)) / 2;
            elseif uTemp(j, 3) > searchRange(2)
                uTemp(j, 3) = (searchRange(2) + x(j, i)) / 2;
            end
        end
        
        % 挑选最优个体
        uTempCost = fhd(uTemp, funcNum);
        [uCost(i), minIndex] = min(uTempCost);
        u(i) = uTemp(:, minIndex);
    end

    % 选择
    goodIndex = uCost <= xCost;
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    trace(g + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);

end

