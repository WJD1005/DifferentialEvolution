clear
clc

mypar = parpool("local", 64);

tic

DE_ConvergenceTest_Parallel_Mat
JADE_ConvergenceTest_Parallel_Mat
CoDE_ConvergenceTest_Parallel_Mat
SHADE_ConvergenceTest_Parallel_Mat
LSHADE_ConvergenceTest_Parallel_Mat
% jSO_ConvergenceTest_Parallel_Mat

toc

delete(mypar);
