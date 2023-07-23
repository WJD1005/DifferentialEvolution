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
xCost = fhd(x, funcNum);  % 初始成本
u = zeros(D, 1);  % 储存单个试验个体

trace = zeros(1, G + 1);  % 储存每代最小值
trace(1) = min(xCost);

% 迭代
for g = 1 : G
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
        if uCost <= xCost(i)
            x(:, i) = u;
            xCost(i) = uCost;
        end
    end
    trace(g + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);
end
