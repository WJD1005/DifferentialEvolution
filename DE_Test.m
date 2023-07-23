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
xCost = fhd(x, funcNum);  % 初始成本
FES = NP;
u = zeros(D, 1);  % 储存单个试验个体

G = ceil(maxFES / NP) - 1;  % 最大代数预分配内存
errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = min(xCost) - realMinVal;

% 迭代
g = 1;
while FES < maxFES
    % DE/rand/1/bin
    for i = 1 : NP
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

        % 保证至少交叉一个维度
        jRand = randi(D);
        
        % 变异交叉
        for j = 1 : D
            if rand() <= CR || j == jRand
                u(j) = x(j, r1) + F * (x(j, r2) - x(j, r3));
                % 越界调整
                if u(j) < searchRange(1)
                    u(j) = (searchRange(1) + x(j, i)) / 2;
                elseif u(j) > searchRange(2)
                    u(j) = (searchRange(2) + x(j, i)) / 2;
                end
            else
                u(j) = x(j, i);
            end
        end

        % 选择
        uCost = fhd(u, funcNum);
        FES = FES + 1;
        if uCost <= xCost(i)
            x(:, i) = u;
            xCost(i) = uCost;
        end

        % 中途强行跳出
        if FES == maxFES
            break
        end
    end
    errorTrace(g + 1) = min(xCost) - realMinVal;
    g = g + 1;
end

minError = errorTrace(end);

end
