% Title: global optimization test
% Author: Anthony Su
% Date: 2024-01-29

close all
clear all
clc

% reproducible rng
rng default

%
gs = GlobalSearch;
objectiveFunction = @(x)(4*x(1)^2 - 2.1*x(1)^4 + x(1)^6/3 ...
    + x(1)*x(2) - 4*x(2)^2 + 4*x(2)^4);
prob = createOptimProblem('fmincon','x0',[-1,2],...
    'objective',objectiveFunction,'lb',[-3,-3],'ub',[3,3]);
[x,fx] = run(gs,prob);