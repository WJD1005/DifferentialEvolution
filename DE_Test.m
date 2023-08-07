function [minError, errorTrace] = DE_Test(NP, D, maxFES, F, CR, searchRange, fhd, funcNum, realMinVal)
%DE_TEST DE算法测试函数，用于测试算法性能。
% 输入：
% NP：种群数量，D：维度，maxFES：最大函数评估次数，F：变异系数，CR：交叉概率，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号，
% realMinVal：真实最小值。
% 输出：
% [minError, errorTrace]
% minError：最小误差值，
% errorTrace：每一代最小误差值记录（1*(G+1))。


% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);
xCost = fhd(x, funcNum);  % 初始成本
uCost = zeros(1, NP);
FES = NP;

G = ceil(maxFES / NP) - 1;  % 最大代数预分配内存
errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = min(xCost) - realMinVal;

% 迭代
g = 1;
while FES < maxFES
    % DE/rand/1/bin
    for i = 1 : NP
        r = randperm(NP, 3);
        regenIndex = r == i;  % 只可能有1个等于i
        while r(regenIndex) == i
            r(regenIndex) = randi(NP);
        end

        % 保证至少交叉一个维度
        jRand = randi(D);
        
        % 变异交叉
        for j = 1 : D
            if rand() <= CR || j == jRand
                u(j, i) = x(j, r(1)) + F * (x(j, r(2)) - x(j, r(3)));
                % 越界调整
                if u(j, i) < searchRange(1)
                    u(j, i) = (searchRange(1) + x(j, i)) / 2;
                elseif u(j, i) > searchRange(2)
                    u(j, i) = (searchRange(2) + x(j, i)) / 2;
                end
            else
                u(j, i) = x(j, i);
            end
        end
    end

    % 选择
    if FES + NP <= maxFES
        uCost = fhd(u, funcNum);
        goodIndex = uCost <= xCost;
        FES = FES + NP;
    else
        uCost(1 : maxFES - FES) = fhd(u(:, 1 : maxFES - FES), funcNum);
        goodIndex = false(1, NP);  % 保证逻辑索引长度
        goodIndex(1 : maxFES - FES) = uCost(1 : maxFES - FES) <= xCost(1 : maxFES - FES);
        FES = maxFES;
    end
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    errorTrace(g + 1) = min(xCost) - realMinVal;
    g = g + 1;
end

minError = errorTrace(end);

end
