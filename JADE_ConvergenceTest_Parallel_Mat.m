clear;
clc;

% 算法参数
NP = 100;
p = 0.05;
c = 0.1;
searchRange = [-100, 100];

% 测试参数
D = [30, 50];  % 测试维度（2/10/30/50/100，可输入向量）
maxFES = 1e4 .* D;  % 最大函数评估次数
fhd = str2func('cec17_func');  % 调用CEC17标准测试集
funcNum = [1, 3 : 30];  % 测试函数序号（可输入向量）
realMinVal = funcNum .* 100;  % 真正最小值
testNum = 51;  % 测试次数

% 保存设置
fileName = 'JADE.mat';  % 文件名

% 测试结果
minError = zeros(length(funcNum), testNum, length(D));  % 最小误差

% 开启并行计算环境
% par = parpool;  % 自动分配核心数

% 多维度测试
for i = 1 : length(D)
    % 多函数测试
    for j = 1 : length(funcNum)
        % 多次测试
        parfor k = 1 : testNum
            [minError(j, k, i), ~] = JADE_Test(NP, D(i), maxFES(i), p, c, searchRange, fhd, funcNum(j), realMinVal(j));
        end
    end
end

% 保存数据
save(fileName, 'minError');

% 关闭并行计算环境
% delete(par);
