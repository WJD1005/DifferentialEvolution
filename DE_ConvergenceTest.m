clear;
clc;

% 算法参数
F = 0.5;
CR = 0.1;
searchRange = [-100, 100];

% 测试参数
D = [2, 10];  % 测试维度
NP = [20, 100];  % 种群数（与测试维度对应）
maxG = 1e5;  % 最大测试代数
fhd = str2func('cec17_func');  % 调用CEC17标准测试集
FuncNum = 1;  % 测试函数序号
testNum = 10;  % 测试次数
realMinVal = 100;  % 真正最小值
errorRange = 1e-6;  % 达到该误差范围即算作结束

% 测试结果
convergenceGen = zeros(length(D), testNum);  % 收敛代数

% 多维度测试
for i = 1 : length(D)
    % 多次测试
    for j = 1 : testNum
        [convergenceGen(i, j), trace] = DE_Test(NP(i), D(i), maxG, F, CR, searchRange, fhd, FuncNum, realMinVal, errorRange);
    end
end


% 无法收敛时画图观察是代数不够还是收敛于局部最小值
% plot(trace)


