% Title: MARGE Tuning Optimization Driver
% Author: Anthony Su
% Date: 2023-10-25

close all
clear all
clc


% ZZ1 = [2,0,1,0,0;
%     2,0,1,1,0;
%     2,0,1,1,1;
%     5,0,1,1,1;
%     2,2,1,1,1;
%     5,2,1,1,1];
% 
% ZZ2 = [2,0,1e3,0,0;
%     2,0,1e3,1,0;
%     2,0,1e3,1,1;
%     5,0,1e3,1,1;
%     2,2,1e3,1,1;
%     5,2,1e3,1,1];
% 
% ZZ3 = [2,0,1e-3,0,0;
%     2,0,1e-3,1,0;
%     2,0,1e-3,1,1;
%     5,0,1e-3,1,1;
%     2,2,1e-3,1,1;
%     5,2,1e-3,1,1];
% 
% for i = 1:6
%     ZZ(3*i-2:3*i,:) = [ZZ1(i,:); ZZ2(i,:); ZZ3(i,:)]
% end

ZZ = [2 2,1,1,0];

for idxOpt = 1:height(ZZ)
    z = ZZ(idxOpt,:);

% initialize 
% mpWeights = logspace(-2,2,9);
ns = z(1)
nLag = z(2)
mpWeights = z(3)
x0_manual = logical(z(4))
infBounds = logical(z(5))

resDatabase = zeros(size(mpWeights));
magErrorDatabase = zeros(size(mpWeights));
phaseErrorDatabase = zeros(size(mpWeights));
xDatabase = zeros(length(mpWeights),1+ns+4+(nLag+3)*ns*ns);

for idxWeight = 1:length(mpWeights)

    %% DEFINE OBJECTIVE FUNCTION
    
    % weight to apply to magnitude error over phase error
    mpWeight = mpWeights(idxWeight);
    
    % anonymous objective function which incorporates the mpWeight
    objectiveFunction = @(x) margeObjective(x,mpWeight,false,ns,nLag);
    
    %% INITIAL CONDITIONS
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
    
    %%%%%%%%%%%%%% +----------------------------------------------------------+
    %% REMEMBER %% | NEED TO REMOVE EXTRAPOLATION BEFORE LOOKING AT ACCELS!!! |
    %%%%%%%%%%%%%% +----------------------------------------------------------+
    
    % inequality constraints v1
    % omega2Lim = 1.454401*[0.9,1.1];
    % zeta1Lim = [0,Inf];
    % zeta2Lim = 0.028*[0.5,1.5];
    % zetaLim = [zeta1Lim',zeta2Lim'] ;
    % dCtrlLim = [0,1];
    % dPLim = [0.5,1.5];

    % inequality constraints v2
    % omega2Lim = 1.454401*[0.9,1.1];
    % zeta1Lim = [0,Inf];
    % zeta2Lim = 0.028*[0.2,5];
    % zetaLim = [zeta1Lim',zeta2Lim'] ;
    % dCtrlLim = [0,5];
    % dPLim = [0,5];

    % inequality constraints v3
    % omega2Lim = 1.454401*[0.9,1.1];
    % zeta1Lim = [0,Inf];
    % zeta2Lim = 0.028*[0.2,5];
    % zetaLim = [zeta1Lim',zeta2Lim'] ;
    % dCtrlLim = [-5,5];
    % dPLim = [-5,5];

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
    
    % % single optimization
    opt = optimoptions('fmincon','UseParallel',true,'Display','final-detailed');
    xNew = fmincon(objectiveFunction,x0,[],[],[],[],LB,UB,[],opt);

    % store result
    xDatabase(idxWeight,:) = xNew';
    [residualNew,magError,phaseError] = objectiveFunction(xNew);
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

%% SAVE MODEL TO ASE_SS.mat

% save final model at equal weighting
% margeObjective(xDatabase(5,:),1,true);

% save final model (when using only one weight)
[residual,magError,phaseError] = margeObjective(xDatabase(1,:),mpWeight,true,ns,nLag);
magError = sum(magError,'all')
phaseError = sum(phaseError,'all')
residual

% save final design vector
bndName = char(infBounds*'Inf'+~infBounds*'Fin');
x0Name = char(x0_manual*'MAN'+~x0_manual*'OG_');
savename = ['optModelParams/ns',num2str(ns),'_nLag',num2str(nLag),'_mpWeight',num2str(log10(mpWeight)),'_bounds',bndName,'_x0',x0Name,'.mat'];
xFinal = xDatabase(1,:);
save(savename,'x0','xFinal','residual','magError','phaseError','residual','ns','nLag','mpWeight','LB','UB')

%% PLOT RESULT

% plot pareto front
% semilogx(mpWeights,magErrorDatabase,'.-','DisplayName','magnitude')
% hold on
% semilogx(mpWeights,phaseErrorDatabase,'.-','DisplayName','phase')
% title('MARGE Optimal Solutions')
% xlabel('ratio of magnitude/phase error weights')
% ylabel('optimal residual')
% grid on
% legend()
% ax = gca;
% ax.YLim(1) = 0;

% plot FRFs comparison
% margeFreqExperiment




XX{idxOpt} = xFinal;

end

% save(['OPTVARS_1e',num2str(log10(mpWeights)),'.mat'],'XX')
save('OPTVARS.mat','XX')