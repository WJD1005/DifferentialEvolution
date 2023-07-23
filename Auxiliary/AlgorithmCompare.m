clear
clc

% B compares with A
A = 'L-SHADE';
B = ["DE", "JADE", "CoDE", "SHADE"];
D = [30, 50];
funcNum = [1, 3 : 30];

fileName = strcat(A, '_Comparison.xls');

result = strings(funcNum(end), 1);

for i = 1 : length(D)
    Adata = xlsread(strcat(A, '.xls'), sprintf('D%d', D(i)), sprintf('A%d:AY%d', funcNum(1), funcNum(end)));
    for j = 1 : length(B)
        Bdata = xlsread(strcat(B(j), '.xls'), sprintf('D%d', D(i)), sprintf('A%d:AY%d', funcNum(1), funcNum(end)));
        for k = 1 : length(funcNum)
            result(funcNum(k)) = ranksumtest(Adata(funcNum(k), :), Bdata(funcNum(k), :));
        end
        xlswrite(fileName, result, sprintf('%s_D%d', B(j), D(i)), 'A1');
    end
end
