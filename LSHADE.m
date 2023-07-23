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
xCost = fhd(x, funcNum);  % 初始成本
FES = NP;
u = zeros(D, 1);  % 储存单个试验个体

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
    
    for i = 1 : NP
        % 生成参数
        ri = randi(H);  % 均值索引
        if mCR(ri) == -1
            CRi = 0;  % 均值为终止值时CR=0
        else
            CRi = normrnd(mCR(ri), 0.1);
        end
        if CRi < 0
            CRi = 0;  % 截断
        elseif CRi > 1
            CRi = 1;  % 截断
        end
        Fi = mF(ri) + 0.1 * tan((rand() - 0.5) * pi);  % 均匀分布转柯西分布
        while Fi <= 0
            Fi = mF(ri) + 0.1 * tan((rand() - 0.5) * pi);  % 重新生成
        end
        if Fi > 1
            Fi = 1;  % 截断
        end
        
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
            if rand() < CRi || j == jRand
                u(j) = x(j, i) + Fi * (x(j, pbest) - x(j, i)) + Fi * (x(j, r1) - PUA(j, r2));
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
            % 避免差为0
            if uCost < xCost(i)
                % 储存失败解
                A(:, end + 1) = x(:, i);
                % 储存成功参数
                SCR(end + 1) = CRi;
                SF(end + 1) = Fi;
                deltaF(end + 1) = xCost(i) - uCost;  % 储存差值
            end

            x(:, i) = u;
            xCost(i) = uCost;
        end

        % 中途强行跳出
        if FES == maxFES
            break
        end
    end
    
    % 保证A的大小不超过Asize
    randomIndex = randperm(size(A, 2));
    A(:, randomIndex(1 : size(A, 2) - Asize)) = [];

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
        xCost(index(1 : NP - nextNP)) = [];  % 删除对应的函数值
        NP = nextNP;
    end

    trace(end + 1) = min(xCost);
end

[minVal, i] = min(xCost);
solution = x(:, i);

end
