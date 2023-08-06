function [minVal, solution, trace] = JADE(NP, D, G, p, c, searchRange, fhd, funcNum)
%JADE JADE算法函数。
% 输入：
% NP：种群数量，D：维度，G：代数，p：优解比例，c：参数自适应权值，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号。
% 输出：
% [minVal, solution, trace]
% minVal：函数最小值，
% solution：函数最优解（D*1）,
% trace：每一代的函数最小值（1*(G+1))。

% 初始参数
uCR = 0.5;
uF = 0.5;
A = [];  % 动态内存方便一点但运行稍慢

% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);
xCost = fhd(x, funcNum);  % 初始成本

trace = zeros(1, G + 1);  % 储存每代最小值
trace(1) = min(xCost);


% 迭代
for g = 1 : G
    SCR = [];
    SF = [];


    % 计算这一代成本并排序
    [~, index] = sort(xCost);

    % P并A
    PUA = [x, A];

    % 生成参数
    CR = uCR + 0.1 .* randn(1, NP);
    CR = min(CR, 1);  % 截断
    CR = max(CR, 0);  % 截断
    F = uF + 0.1 .* tan((rand(1, NP) - 0.5) .* pi);
    regenIndex = find(F <= 0);
    while ~isempty(regenIndex)
        F(regenIndex) = uF + 0.1 .* tan((rand(1, length(regenIndex)) - 0.5) .* pi);  % 重新生成
        regenIndex = find(F <= 0);
    end
    F = min(F, 1);  % 截断
    
    for i = 1 : NP
        % 取优解
        pbest = index(randi(max(round(p * NP), 2)));  % 最少取两个
        
        % 取不重复随机解
        r1 = randi(NP);
        while r1 == i
            r1 = randi(NP);
        end
        r2 = randi(size(PUA, 2));
        while r2 == i || r2 == r1
            r2 = randi(size(PUA, 2));
        end
        
        % 保证至少交叉一个维度
        jRand = randi(D);
        
        % 变异交叉
        for j = 1 : D
            if rand() <= CR(i) || j == jRand
                u(j, i) = x(j, i) + F(i) * (x(j, pbest) - x(j, i)) + F(i) * (x(j, r1) - PUA(j, r2));
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
    % 储存失败解
    A = [A, x(:, goodIndex)];
    % 保证A的大小不超过NP
    if size(A, 2) > NP
        randomIndex = randperm(size(A, 2), NP);
        A = A(:, randomIndex);
    end
    % 储存成功参数
    SCR = CR(goodIndex);
    SF = F(goodIndex);
    % 取代
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    % 集合非空时更新自适应参数均值
    if ~isempty(SCR)
        uCR = (1 - c) * uCR + c * mean(SCR);
        uF = (1 - c) * uF + c * (sum(SF .^ 2) / sum(SF));
    end

    trace(g + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);

end
