function [minVal, solution, trace] = LSHADE(initNP, D, maxFES, p, Asize, H, searchRange, fhd, funcNum)
%LSHADE L-SHADE算法函数。
% 输入：
% initNP：初始种群数量，D：维度，maxFES：最大函数评估次数，
% p：优解比例，Asize：次优解集大小，H：参数均值集大小，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号。
% 输出：
% [minVal, solution, trace]
% minVal：函数最小值，
% solution：函数最优解（D*1）,
% trace：每一代的函数最小值（1*(G+1))。


% 初始参数
NP = initNP;
minNP = 4;  % 种群数量最小值
mCR = ones(1, H) .* 0.5;
mF = ones(1, H) .* 0.5;
k = 1;
A = [];  % 动态内存方便一点但运行稍慢

% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);
xCost = fhd(x, funcNum);  % 初始成本
FES = NP;

trace(1) = min(xCost);  % 储存每代最小值


% 迭代
while FES < maxFES
    SCR = [];
    SF = [];
    deltaF = [];  % 用于储存函数差值便于计算权重

    % 计算这一代成本并排序
    [~, index] = sort(xCost);

    % P并A
    PUA = [x, A];

    % 生成参数
    r = randi(H, [1, NP]);
    CR = mCR(r) + 0.1 .* randn(1, NP);
    CR(mCR(r) == -1) = 0;  % 均值为终止值时CR=0
    CR = min(CR, 1);  % 截断
    CR = max(CR, 0);  % 截断
    F = mF(r) + 0.1 .* tan((rand(1, NP) - 0.5) .* pi);
    regenIndex = find(F <= 0);
    while ~isempty(regenIndex)
        F(regenIndex) = mF(r(regenIndex)) + 0.1 .* tan((rand(1, length(regenIndex)) - 0.5) .* pi);  % 重新生成
        regenIndex = find(F <= 0);
    end
    F = min(F, 1);  % 截断
    
    for i = 1 : NP
        % 取优解
        pbest = index(randi(max(round(p * NP), 2)));  % 最少从前2个里取
        
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
        realGoodIndex = uCost < xCost;
        FES = FES + NP;
    else
        uCost(1 : maxFES - FES) = fhd(u(:, 1 : maxFES - FES), funcNum);
        goodIndex = false(1, NP);  % 保证逻辑索引长度
        realGoodIndex = false(1, NP);  % 保证逻辑索引长度
        goodIndex(1 : maxFES - FES) = uCost(1 : maxFES - FES) <= xCost(1 : maxFES - FES);
        realGoodIndex(1 : maxFES - FES) = uCost(1 : maxFES - FES) < xCost(1 : maxFES - FES);
        FES = maxFES;
    end
    % 储存失败解
    A = [A, x(:, realGoodIndex)];
    % 保证A的大小不超过Asize
    if size(A, 2) > Asize
        randomIndex = randperm(size(A, 2), Asize);
        A = A(:, randomIndex);
    end
    % 储存成功参数
    SCR = CR(realGoodIndex);
    SF = F(realGoodIndex);
    deltaF = xCost(realGoodIndex) - uCost(realGoodIndex);
    % 取代
    x(:, goodIndex) = u(:, goodIndex);
    xCost(goodIndex) = uCost(goodIndex);

    % 集合非空时更新自适应参数均值
    if ~isempty(SCR)
        w = deltaF ./ sum(deltaF);  % 权重
        if mCR(k) == -1 || max(SCR) == 0
            mCR(k) = -1;  % 置终止值
        else
            mCR(k) = sum(w .* SCR .^ 2) / sum(w .* SCR);
        end
        mF(k) = sum(w .* SF .^ 2) / sum(w .* SF);

        k = k + 1;
        if k > H
            k = 1;
        end
    end
    
    % LPSR
    nextNP = round((minNP - initNP) / maxFES * FES + initNP);  % 下一代种群数量
    if nextNP < NP
        [~, index] = sort(xCost, 'descend');  % 降序排序索引
        x(:, index(1 : NP - nextNP)) = [];  % 删除个体
        u(:, index(1 : NP - nextNP)) = [];
        xCost(index(1 : NP - nextNP)) = [];  % 删除对应的函数值
        uCost(index(1 : NP - nextNP)) = [];
        NP = nextNP;
    end

    trace(end + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);

end
