function [convergenceGen, trace] = JADE_Test(NP, D, maxG, p, c, searchRange, fhd, funcNum, realMinVal, errorRange)
%JADE JADE算法测试函数，用于测试收敛性。
% 输入：
% NP：种群数量，D：维度，maxG：最大代数，p：优解比例，c：参数自适应权值，
% searchRange：搜索范围（1*2），fhd：测试函数句柄，funcNum：测试函数序号，
% realMinVal：真实最小值，errorRange：误差范围。
% 输出：
% [convergenceGen, trace]
% convergenceGen：最小值与真实值误差在误差范围内时的代数，
% trace：每一代的函数最小值（1*(convergenceGen+1))。

% 初始参数
uCR = 0.5;
uF = 0.5;
A = [];  % 动态内存方便一点但运行稍慢

% 初始种群
x = rand(D, NP) .* (searchRange(2) - searchRange(1)) + searchRange(1);
u = zeros(D, NP);

trace = zeros(1, maxG + 1);  % 储存每代最小值
trace(1) = min(fhd(x, funcNum));
convergenceGen = 0;


% 迭代
for g = 1 : maxG
    SCR = [];
    SF = [];

    % 计算这一代成本并排序
    xCost = fhd(x, funcNum);
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
        pbest = index(randi(ceil(p * NP)));  % 采用向上取整
        
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
                u(j, i) = x(j, i) + Fi * (x(j, pbest) - x(j, i)) + Fi * (x(j, r1) - PUA(j, r2));
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
            % 储存失败解
            A(:, end + 1) = x(:, i);

            % 储存成功参数
            SCR(end + 1) = CRi;
            SF(end + 1) = Fi;

            x(:, i) = u(:, i);
        end
    end
    % 保证A的大小不超过NP
    while size(A, 2) > NP
        A(:, randi(size(A, 2))) = [];
    end

    % 更新自适应参数均值
    uCR = (1 - c) * uCR + c * mean(SCR);
    uF = (1 - c) * uF + c * (sum(SF .^ 2) / sum(SF));

    trace(g + 1) = min(fhd(x, funcNum));

    % 判断是否在误差范围内
    if abs(trace(g + 1) - realMinVal) <= errorRange
        trace = trace(1 : g + 1);  % 截断
        convergenceGen = g;
        break
    end
end

end
