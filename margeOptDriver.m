% Title: MARGE Tuning Optimization Driver
% Author: Anthony Su
% Date: 2023-10-25

close all
clear all
clc

% initialize 
mpWeights = logspace(-1,3,9);
resDatabase = zeros(size(mpWeights));
magErrorDatabase = zeros(size(mpWeights));
phaseErrorDatabase = zeros(size(mpWeights));
xDatabase = zeros(length(mpWeights),19);

for idxWeight = 1:length(mpWeights)

    %% DEFINE OBJECTIVE FUNCTION
    
    % weight to apply to magnitude error over phase error
    mpWeight = mpWeights(idxWeight);
    
    % anonymous objective function which incorporates the mpWeight
    objectiveFunction = @(x) margeObjective(x,mpWeight);
    
    %% INITIAL CONDITIONS
    
    % initialize design variables
    omega2 = 1.454401; % Hz, set to zero to ignore
    zeta1 = 0.3*200;
    zeta2 = 0.028*1;
    dAil1 = 0.6;
    dAil2 = 0.7;
    dElev = 0.6;
    % dVane = 4;
    dVane = 1;
    dP = ones(2,2,3);
    dP(1,1,:) = 0.9;
    dP(1,2,:) = 1;
    dP(2,1,:) = 0.5;
    dP(2,2,:) = 1.5;
    
    x0 = [omega2,zeta1,zeta2,dAil1,dAil2,dElev,dVane,reshape(dP,1,[])];
    
    % compute initial residual
    [residual0,magError0,phaseError0] = objectiveFunction(x0);
    
    %%%%%%%%%%%%%% +----------------------------------------------------------+
    %% REMEMBER %% | NEED TO REMOVE EXTRAPOLATION BEFORE LOOKING AT ACCELS!!! |
    %%%%%%%%%%%%%% +----------------------------------------------------------+
    
    % inequality constraints
    omega2Lim = 1.454401*[0.9,1.1];
    zeta1Lim = [0,Inf];
    zeta2Lim = 0.028*[0.5,1.5];
    dCtrlLim = [0,1];
    dPLim = [0.5,1.5];
    
    LB = [omega2Lim(1),zeta1Lim(1),zeta2Lim(1),dCtrlLim(1).*ones(1,4),dPLim(1).*ones(1,2*2*3)];
    UB = [omega2Lim(2),zeta1Lim(2),zeta2Lim(2),dCtrlLim(2).*ones(1,4),dPLim(2).*ones(1,2*2*3)];
    
    % optimization
    opt = optimoptions('fmincon','UseParallel',true,'Display','final-detailed');
    xNew = fmincon(objectiveFunction,x0,[],[],[],[],LB,UB,[],opt);

    % store result
    xDatabase(idxWeight,:) = xNew';
    [residualNew,magError,phaseError] = margeObjective(xNew,mpWeight);
    magErrorDatabase(idxWeight) = norm(magError,'fro');
    phaseErrorDatabase(idxWeight) = norm(phaseError,'fro');

    % print result
    % varNames = {'omega2','zeta1','zeta2','dAil1','dAil2','dElev','dVane',...
    %     'dP111','dP211','dP121','dP221','dP112','dP212','dP122','dP222',...
    %     'dP113','dP213','dP123','dP223'};
    % fprintf('improved from initial residual %6.3f to final residual %6.3f\n',residual0,residualNew)
    % fprintf('   varName     x0   xNew\n')
    % for idxVar = 1:length(x0)
    % fprintf('%10s %6.3f %6.3f\n',varNames{idxVar},x0(idxVar),xNew(idxVar))
    % end

end

%% PLOT RESULT

% plot pareto front
semilogx(mpWeights,magErrorDatabase,'.-')
hold on
semilogx(mpWeights,phaseErrorDatabase,'.-')
title('MARGE Optimal Solutions')
xlabel('ratio of magnitude/phase error weights')
ylabel('optimal residual')
grid on