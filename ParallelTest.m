clear
clc

mypar = parpool("local", 64);

tic

DE_ConvergenceTest_Parallel
JADE_ConvergenceTest_Parallel
CoDE_ConvergenceTest_Parallel
SHADE_ConvergenceTest_Parallel
LSHADE_ConvergenceTest_Parallel
% jSO_ConvergenceTest_Parallel

toc

delete(mypar);
