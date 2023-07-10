function [convergenceGen, trace] = DE_Test(NP, D, maxG, F, CR, searchRange, fhd, funcNum, realMinVal, errorRange)
%DE DE算法测试函数，用于测试收敛性。
% 输入：
% NP：种群数量，D：维度，maxG：最大代数，F：变异系数，CR：交叉系数，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号，
% realMinVal：真实最小值，errorRange：误差范围。
% 输出：
% [convergenceGen, trace]
% convergenceGen：最小值与真实值误差在误差范围内时的代数，
% trace：每一代的函数最小值（1*(G+1))。


% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);
trace = zeros(1, maxG + 1);  % 储存每代最小值
trace(1) = min(fhd(x, funcNum));
convergenceGen = 0;

% 迭代
for g = 1 : maxG
    % DE/rand/1/bin
    for i = 1:NP
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

        for j = 1 : D
            if rand() <= CR || j == jRand
                u(j, i) = x(j, r1) + F * (x(j, r2) - x(j, r3));
            else
                u(j, i) = x(j, i);
            end

            % 越界截断
            if u(j, i) < searchRange(1)
                u(j, i) = searchRange(1);
            elseif u(j, i) > searchRange(2)
                u(j, i) = searchRange(2);
            end
        end
    end

    % 选择
    xCost = fhd(x, funcNum);
    uCost = fhd(u, funcNum);
    for i = 1 : NP
        if uCost(i) <= xCost(i)
            x(:, i) = u(:, i);
        end
    end
    trace(g + 1) = min(fhd(x, funcNum));

    % 判断是否在误差范围内
    if abs(trace(g + 1) - realMinVal) <= errorRange
        trace = trace(1, g + 1);  % 截断
        convergenceGen = g;
        break
    end
end

end
