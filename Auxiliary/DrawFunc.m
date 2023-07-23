clear
clc
close all

fhd = str2func('cec17_func');
funcNum = 6;

x = -100 : 100;
y = -100 : 100;
[X, Y] = meshgrid(x, y);
x = repelem(x, 201);
y = repmat(y, 1, 201);
xy = [x; y];
Z = fhd(xy, funcNum);
Z = reshape(Z, [201, 201]);
surf(X,Y,Z);



