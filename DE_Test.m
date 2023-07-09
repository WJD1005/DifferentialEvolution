clear;
clc;

% 算法参数
NP = 50;
D = 10;
G = 100000;
F = 0.5;
CR = 0.1;
searchRange = [-100, 100];

% 测试参数
fhd = str2func('cec17_func');  % 调用CEC17标准测试集
FuncNum = 1;  % 测试函数序号
testNum = 10;  % 测试次数
realMinVal = 100;  % 真正最小值
errorRange = 1e-06;  % 达到该误差范围即算作结束

% 测试结果
iterations = zeros(1, testNum);  % 迭代次数


for i = 1 : testNum
    [minVal, solution, trace] = DE(NP, D, G, F, CR, searchRange, fhd, FuncNum);

    % 记录达到误差范围的迭代次数
    error = trace - realMinVal;
    index = find(abs(error) <= errorRange, 1, "first");
    iterations(i) = index - 1;
end


