% Title: MARGE Tuning Global Optimization Driver
% Author: Anthony Su
% Date: 2024-01-28

close all
clear all
clc

rng default % enforce repeatability

% define optimization/model meta-parameters
ns = 2
nLag = 2
mpWeight = 1 % magnitude weight vs phase weight
x0_manual = true % whether to use manually tuned initial guess (as opposed to un-tuned)
infBounds = false % whether to not bound solution to physical values

%% DEFINE OBJECTIVE FUNCTION

% anonymous objective function which incorporates the mpWeight
objectiveFunction = @(x) margeObjective(x,mpWeight,false,ns,nLag);

%% OPTIMIZATION INITIAL CONDITIONS
zeta(1:9) = [0,0.028,0.042,0.030,0.112,0.031,0.018,0.022,0.032];

if(~x0_manual)
% initialize design variables to untuned result
    omega2 = 1.454401; % Hz, set to zero to ignore
    zeta(1) = 0;
    zeta(2) = 0.028*1;
    dAil1 = 1;
    dAil2 = 1;
    dElev = 1;
    dVane = 1;
    dP = ones(ns,ns,3+nLag);
    dP(1,1,:) = 1;
    dP(1,2,:) = 1;
    dP(2,1,:) = 1;
    dP(2,2,:) = 1;
else
    % initialize design variables to maually tuned result
    omega2 = 1.454401; % Hz, set to zero to ignore
    zeta(1) = 0.3*200;
    zeta(2) = 0.028*1;
    dAil1 = 0.6;
    dAil2 = 0.7;
    dElev = 0.6;
    % dVane = 4;
    dVane = 1;
    dP = ones(ns,ns,3+nLag);
    dP(1,1,:) = 0.9;
    dP(1,2,:) = 1;
    dP(2,1,:) = 0.5;
    dP(2,2,:) = 1.5;
end

zeta = zeta(1:ns);
x0 = [omega2,zeta,dAil1,dAil2,dElev,dVane,reshape(dP,1,[])];

% compute initial residual
[residual0,magError0,phaseError0] = objectiveFunction(x0);

%% OPTIMIZATION BOUNDS

if(~infBounds)
    % inequality constraints v4 (used in optimization studies)
    omega2Lim = 1.454401*[0.9,1.1];
    zeta1Lim = [0,Inf];
    zeta2Lim = 0.028*[0.5,1.5];
    zetaLim = [zeta1Lim',zeta2Lim'] ;
    dCtrlLim = [0,1];
    dPLim = [0,Inf];
else
    % no inequality constraints (used in optimization studies)
    omega2Lim = [-Inf,Inf];
    zetaLim = [-Inf;Inf].*ones(1,ns);
    dCtrlLim = [-Inf,Inf];
    dPLim = [-Inf,Inf];
end

LB = [omega2Lim(1),zetaLim(1,:),dCtrlLim(1).*ones(1,4),dPLim(1).*ones(1,ns*ns*(3+nLag))];
UB = [omega2Lim(2),zetaLim(2,:),dCtrlLim(2).*ones(1,4),dPLim(2).*ones(1,ns*ns*(3+nLag))];

%% OPTIMIZE

% % single optimization
% opt = optimoptions('fmincon','UseParallel',true,'Display','final-detailed');
% xNew = fmincon(objectiveFunction,x0,[],[],[],[],LB,UB,[],opt);

% % multiple optimization ("global")
% xVar = optimvar("x",length(x0),LowerBound=LB,UpperBound=UB);
% x0_.xVar = x0;
% prob = optimproblem(Objective=objectiveFunction(xVar));
% opt = optimdoptions('fmincon','UseParallel',true,'Display','final-detailed');
% ms = GlobalSearch;
% [xNew,fval] = solve(prob,x0_,options=opt);

% multiple optimization ("global") try #2
gs = GlobalSearch('Display','iter','MaxTime',5000); % ~1.4 hour runtime
% gs = GlobalSearch('Display','iter','MaxTime',1000); % 16min 40sec hour runtime
opt = optimoptions('fmincon',UseParallel=true);
problem = createOptimProblem('fmincon',x0=x0,...
objective=objectiveFunction,lb=LB,ub=UB,options=opt);
xFinal = run(gs,problem);

%%

% store result
[residual,magError,phaseError] = objectiveFunction(xFinal);
magError = sum(magError,'all');
phaseError = sum(phaseError,'all');

%% SAVE MODEL TO ASE_SS.mat

wantSave = input("Save result? (Y/N): ",'s');
if(wantSave=='Y')
    wantSave = true;
    [residual,magError,phaseError] = margeObjective(xFinal,mpWeight,true,ns,nLag);
elseif(wantSave=='N')
    wantSave = false;
else
    error('wantSave must be Y or N')
end

% save final design vector
if(wantSave)
    bndName = char(infBounds*'Inf'+~infBounds*'Fin');
    x0Name = char(x0_manual*'MAN'+~x0_manual*'OG_');
    savename = ['optModelParams/ns',num2str(ns),'_nLag',num2str(nLag),'_mpWeight',num2str(log10(mpWeight)),'_bounds',bndName,'_x0',x0Name,'_GLOBAL.mat'];
    save(savename,'x0','xFinal','residual','magError','phaseError','ns','nLag','mpWeight','LB','UB')
end

%% PLOT RESULT

% plot FRFs comparison
margeFreqExperiment
