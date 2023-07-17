clear;
clc;

% 算法参数
NP = 100;
p = 0.05;
c = 0.1;
searchRange = [-100, 100];

% 测试参数
D = 30;  % 测试维度（2/10/30/50/100）
G = 1e4 * D / NP;  % 测试代数
fhd = str2func('cec17_func');  % 调用CEC17标准测试集
funcNum = [1, 3 : 30];  % 测试函数序号（可输入向量）
realMinVal = funcNum .* 100;  % 真正最小值
testNum = 51;  % 测试次数

% 保存设置
fileName = 'JADE.xls';  % 文件名
rawDataSheetName = sprintf('D%d', D);  % 原始数据工作表名
statisticalDataSheetName = sprintf('D%dmean&std', D);  % 统计数据工作表名

% 单函数测试结果
minError = zeros(1, testNum);  % 最小误差

% 多函数测试
for i = 1 : length(funcNum)
    % 多次测试
    for j = 1 : testNum
        [minError(j), ~] = JADE_Test(NP, D, G, p, c, searchRange, fhd, funcNum(i), realMinVal(i));
    end

    % 保存结果
    origin = sprintf('A%d', funcNum(i));  % 数据起始单元格
    xlswrite(fileName, sort(minError), rawDataSheetName, origin);  % 原始数据
    xlswrite(fileName, [mean(minError), std(minError)], statisticalDataSheetName, origin);  % 统计数据
end
