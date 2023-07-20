function [minError, errorTrace] = JADE_Test(NP, D, maxFES, p, c, searchRange, fhd, funcNum, realMinVal)
%JADE JADE算法测试函数，用于测试算法性能。
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
xCost = fhd(x, funcNum);  % 初始成本
FES = NP;
u = zeros(D, 1);  % 储存单个试验个体

G = ceil(maxFES / NP) - 1;  % 最大代数预分配内存
errorTrace = zeros(1, G + 1);  % 储存每代最小误差值
errorTrace(1) = abs(min(xCost) - realMinVal);


% 迭代
g = 1;
while FES < maxFES
    SCR = [];
    SF = [];

    % 计算这一代成本并排序
    [~, index] = sort(xCost);

    % P并A
    PUA = [x, A];
    
    for i = 1 : NP
        % 生成参数
        CRi = normrnd(uCR, 0.1);
        if CRi < 0
            CRi = 0;  % 截断
        elseif CRi > 1
            CRi = 1;  % 截断
        end

        Fi = uF + 0.1 * tan((rand() - 0.5) * pi);  % 均匀分布转柯西分布
        while Fi <= 0
            Fi = uF + 0.1 * tan((rand() - 0.5) * pi);  % 重新生成
        end
        if Fi > 1
            Fi = 1;  % 截断
        end

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
            if rand() < CRi || j == jRand
                u(j) = x(j, i) + Fi * (x(j, pbest) - x(j, i)) + Fi * (x(j, r1) - PUA(j, r2));
            else
                u(j) = x(j, i);
            end

            % 越界截断
            if u(j) < searchRange(1)
                u(j) = searchRange(1);
            elseif u(j) > searchRange(2)
                u(j) = searchRange(2);
            end
        end
        
        % 选择
        uCost = fhd(u, funcNum);
        FES = FES + 1;
        if uCost <= xCost(i)
            % 储存失败解
            A(:, end + 1) = x(:, i);

            % 储存成功参数
            SCR(end + 1) = CRi;
            SF(end + 1) = Fi;

            x(:, i) = u;
            xCost(i) = uCost;
        end

        % 中途强行跳出
        if FES == maxFES
            break
        end
    end
    % 保证A的大小不超过NP
    while size(A, 2) > NP
        A(:, randi(size(A, 2))) = [];
    end

    % 集合非空时更新自适应参数均值
    if ~isempty(SCR)
        uCR = (1 - c) * uCR + c * mean(SCR);
        uF = (1 - c) * uF + c * (sum(SF .^ 2) / sum(SF));
    end

    errorTrace(g + 1) = abs(min(xCost) - realMinVal);
    g = g + 1;
end

minError = errorTrace(end);

end
