clear
clc

% 使用作者提供代码

% 测试参数
D = [30, 50];  % 测试维度（2/10/30/50/100，可输入向量）
funcNum = [1, 3 : 30];  % 测试函数序号（可输入向量）
testNum = 51;  % 测试次数

% 保存设置
fileName = 'jSO.mat';  % 文件名

% 测试结果
minError = zeros(length(funcNum), testNum, length(D));  % 最小误差

% 多维度测试
for i = 1 : length(D)
    % 多函数测试
    for j = 1 : length(funcNum)
        minError(j, :, i) = jSO(funcNum(j), testNum, D(i));
    end
end

