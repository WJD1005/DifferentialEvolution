function [minVal, solution, trace] = DE(NP, D, G, F, CR, searchRange, fhd, funcNum)
%DE DE算法函数。
% 输入：
% NP：种群数量，D：维度，G：代数，F：变异系数，CR：交叉系数，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号。
% 输出：
% [minVal, solution, trace]
% minVal：函数最小值，
% solution：函数最优解（D*1）,
% trace：每一代的函数最小值（1*(G+1))。


% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);
xCost = fhd(x, funcNum);  % 初始成本

trace = zeros(1, G + 1);  % 储存每代最小值
trace(1) = min(xCost);

% 迭代
for g = 1 : G
    % DE/rand/1/bin
    for i = 1 : NP
        r = randperm(NP, 3);
        regenIndex = find(r == i);
        while ~isempty(regenIndex)
            r(regenIndex) = randi(NP, [1, length(regenIndex)]);
            regenIndex = find(r == i);
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
    uCost = fhd(u, funcNum);
    goodIndex = uCost <= xCost;
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    trace(g + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);
end
