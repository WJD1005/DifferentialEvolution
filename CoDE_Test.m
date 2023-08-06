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
u = zeros(D, NP);
uTemp = zeros(D, 3);  % 暂存三个试验个体
xCost = fhd(x, funcNum);  % 初始成本
uCost = zeros(1, NP);
endi = -1;
FES = NP;

G = ceil((maxFES - NP) / NP / 3);  % 最大代数预分配内存
errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = min(xCost) - realMinVal;

% 迭代
g = 1;
while FES < maxFES
    for i = 1 : NP
        % DE/rand/1/bin
        parameterIndex = randi(3);  % 选择参数组合
        r = randperm(NP, 3);
        regenIndex = find(r == i);
        while ~isempty(regenIndex)
            r(regenIndex) = randi(NP, [1, length(regenIndex)]);
            regenIndex = find(r == i);
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
        regenIndex = find(r == i);
        while ~isempty(regenIndex)
            r(regenIndex) = randi(NP, [1, length(regenIndex)]);
            regenIndex = find(r == i);
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
        regenIndex = find(r == i);
        while ~isempty(regenIndex)
            r(regenIndex) = randi(NP, [1, length(regenIndex)]);
            regenIndex = find(r == i);
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
        if FES + 3 < maxFES
            uTempCost = fhd(uTemp, funcNum);
            [uCost(i), minIndex] = min(uTempCost);
            u(:, i) = uTemp(:, minIndex);
            FES = FES + 3;
        else
            uTempCost = fhd(uTemp(:, 1 : maxFES - FES), funcNum);
            [uCost(i), minIndex] = min(uTempCost);
            u(:, i) = uTemp(:, minIndex);
            endi = i;  % 标记中止个体
            FES = maxFES;
            break;
        end
    end

    % 选择
    if endi == -1
        goodIndex = uCost <= xCost;
    else
        goodIndex = false(1, NP);  % 保证逻辑索引长度
        goodIndex(1 : endi) = uCost(1 : endi) <= xCost(1 : endi);
    end
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    errorTrace(g + 1) = min(xCost) - realMinVal;
    g = g + 1;
end

minError = errorTrace(end);

end
