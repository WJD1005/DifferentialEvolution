clear;
clc;

% 算法参数
F = 0.5;
CR = 0.1;
searchRange = [-100, 100];

% 测试参数
D = 10;  % 测试维度（2/10/20/30/50/100）
NP = 7.5 * D;  % 种群数
maxG = 1e5;  % 最大测试代数
fhd = str2func('cec17_func');  % 调用CEC17标准测试集
funcNum = 30;  % 测试函数序号（可输入向量）
realMinVal = funcNum .* 100;  % 真正最小值
errorRange = 1e-6;  % 达到该误差范围即算作结束
testNum = 5;  % 测试次数

% 测试结果
convergenceGen = zeros(length(funcNum), testNum);  % 收敛代数

% 多函数测试
for i = 1 : length(funcNum)
    % 多次测试
    for j = 1 : testNum
        [convergenceGen(i, j), trace] = DE_Test(NP, D, maxG, F, CR, searchRange, fhd, funcNum(i), realMinVal(i), errorRange);
    end
end

% 保存结果
save('Result.mat', 'convergenceGen');

% 无法收敛时画图观察是代数不够还是收敛于局部最小值
% plot(trace)


