function [minError, errorTrace] = JADE_Test(NP, D, maxFES, p, c, searchRange, fhd, funcNum, realMinVal)
%JADE_TEST JADE算法测试函数，用于测试算法性能。
% 输入：
% NP：种群数量，D：维度，maxFES：最大函数评估次数，p：优解比例，c：参数自适应权值，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号，
% realMinVal：真实最小值。
% 输出：
% [minError, errorTrace]
% minError：最小误差值，
% errorTrace：每一代最小误差值记录（1*(G+1))。

% 初始参数
uCR = 0.5;
uF = 0.5;
A = [];  % 动态内存方便一点但运行稍慢

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

    errorTrace(g + 1) = min(xCost) - realMinVal;
    g = g + 1;
end

minError = errorTrace(end);

end
