clear
clc

input = ["DE", "JADE", "CoDE", "SHADE", "L-SHADE", "jSO"];
D = [30, 50];
funcNum = [1, 3 : 30];

% 多文件
for i = 1 : length(input)
    output = strcat(input(i), '.xls');  % 输出文件名
    load(strcat(input(i), '.mat'));  % 加载原始数据
    % 多维度
    for j = 1 : length(D)
        % 定义维度工作表
        rawDataSheetName = sprintf('D%d', D(j));  % 原始数据工作表名
        statisticalDataSheetName = sprintf('D%dmean&std', D(j));  % 统计数据工作表名
        % 多函数
        for k = 1 : length(funcNum)
            origin = sprintf('A%d', funcNum(k));  % 数据起始单元格
            xlswrite(output, sort(minError(k, :, j)), rawDataSheetName, origin);  % 原始数据
            xlswrite(output, [mean(minError(k, :, j)), std(minError(k, :, j))], statisticalDataSheetName, origin);  % 统计数据
        end
    end
end
