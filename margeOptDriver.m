% Title: MARGE Tuning Optimization Driver
% Author: Anthony Su
% Date: 2023-10-25

close all
clear all
clc

% initialize design variables
omega2 = 1.454; % Hz
zeta1 = 0.3*200;
zeta2 = 0.028*1;
dAil1 = 0.6;
dAil2 = 0.7;
dElev = 0.6;
dVane = 4;
dP = ones(2,2,3);
dP(1,1,:) = 0.9;
dP(1,2,:) = 1;
dP(2,1,:) = 0.5;
dP(2,2,:) = 1.5;

x0 = [omega2,zeta1,zeta2,dAil1,dAil2,dElev,dVane,reshape(dP,1,[])];

% compute initial residual
residual = margeObjective(x0)