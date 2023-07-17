function [minError, errorTrace] = DE_Test(NP, D, G, F, CR, searchRange, fhd, funcNum, realMinVal)
%DE_Test DE算法测试函数，用于测试算法性能。
% 输入：
% NP：种群数量，D：维度，G：代数，F：变异系数，CR：交叉概率，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号，
% realMinVal：真实最小值。
% 输出：
% [minError, errorTrace]
% minError：最小误差值，
% errorTrace：每一代最小误差值记录（1*(G+1))。


% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);

errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = abs(min(fhd(x, funcNum)) - realMinVal);

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

        % 选择
        xCost = fhd(x(:, i), funcNum);
        uCost = fhd(u(:, i), funcNum);
        if uCost <= xCost
            x(:, i) = u(:, i);
        end
    end
    errorTrace(g + 1) = abs(min(fhd(x, funcNum)) - realMinVal);
end

minError = errorTrace(end);

end
