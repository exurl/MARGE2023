% Title: Compare MARGE Time Response to Experimental Data
% Author: Anthony Su
% Date: 2023-09-16

clear all
close all

%% IMPORT STATE-SPACE MODEL

aseModel = load('./ASE_SS.mat');
A = aseModel.A;
Bc = aseModel.Bc;
C = aseModel.C;
Dc = aseModel.Dc;

% for full-size model:
    % nStates
        % 15 flexible modes
        % 15 flexible mode derivatives
        % 15 flexible mode aero lags (1)
        % 15 flexible mode aero lags (2)
        % 15 flexible mode aero lags (3)
        % 15 flexible mode aero lags (4)
        % ~~gust vane aero lag (1-32)~~
    % nInputs
        % 4 control inputs [ail1, ail2, elev, gust]
    % nOutputs
        % 3 accelerometer outputs [2223,2423,3008]
        % 1 root strain output
        % 1 pitch output
        % 1 pitch output derivative

omegaMax = aseModel.omegaMax;
rho = aseModel.rho;
q = aseModel.q;

% u will represent cotrol input below
clear u

%% REPLICATE EXPERIMENTAL RESPONSE WITH STATE-SPACE MODEL

% q60
expObjs(1,1) = load('windTunnel\wtData\FRF_q60_A1S5.mat').data;
expObjs(2,1) = load('windTunnel\wtData\FRF_q60_A2S5.mat').data;
expObjs(3,1) = load('windTunnel\wtData\FRF_q60_ES2.mat').data;
expObjs(4,1) = load('windTunnel\wtData\FRF_q60_GD4.mat').data;

% q100
expObjs(1,2) = load('windTunnel\wtData\FRF_q100_A1S5.mat').data;
expObjs(2,2) = load('windTunnel\wtData\FRF_q100_A2S5.mat').data;
expObjs(3,2) = load('windTunnel\wtData\FRF_q100_ES2.mat').data;
expObjs(4,2) = load('windTunnel\wtData\FRF_q100_GD4.mat').data;

% q164
expObjs(1,3) = load('windTunnel\wtData\FRF_q164_A1S5.mat').data;
expObjs(2,3) = load('windTunnel\wtData\FRF_q164_A2S5.mat').data;
expObjs(3,3) = load('windTunnel\wtData\FRF_q164_ES2.mat').data;
expObjs(4,3) = load('windTunnel\wtData\FRF_q164_GD4.mat').data;

% q207
expObjs(1,4) = load('windTunnel\wtData\FRF_q207_A1S5.mat').data;
expObjs(2,4) = load('windTunnel\wtData\FRF_q207_A2S5.mat').data;
expObjs(3,4) = load('windTunnel\wtData\FRF_q207_ES2.mat').data;
expObjs(4,4) = load('windTunnel\wtData\FRF_q207_GD4.mat').data;

% q281
expObjs(1,5) = load('windTunnel\wtData\FRF_q281_A1S5.mat').data;
expObjs(2,5) = load('windTunnel\wtData\FRF_q281_A2S5.mat').data;
expObjs(3,5) = load('windTunnel\wtData\FRF_q281_ES2.mat').data;
expObjs(4,5) = load('windTunnel\wtData\FRF_q281_GD4.mat').data;

% q343
expObjs(1,6) = load('windTunnel\wtData\FRF_q343_A1S3p5.mat').data;
expObjs(2,6) = load('windTunnel\wtData\FRF_q343_A2S3p5.mat').data;
expObjs(3,6) = load('windTunnel\wtData\FRF_q343_ES1.mat').data;
expObjs(4,6) = load('windTunnel\wtData\FRF_q343_GD4.mat').data;

% state-space response
tic
waitMessage = parfor_wait(24,'Waitbar',true); % using package "waitbar for parfor" by Yun Pu
parfor idxSpeed = 1:6
    for idxInput = 1:4
        waitMessage.Send;
        expObj = expObjs(idxInput,idxSpeed);
        ssObjs(idxInput,idxSpeed) = ssResponse(expObj,A(:,:,idxSpeed),Bc(:,:,idxSpeed),C(:,:,idxSpeed),Dc(:,:,idxSpeed));
    end
end
waitMessage.Destroy;
toc

%% SAVE EVERYTHING TOGETHER
save('TIME_compare.mat','expObjs','ssObjs','q','rho','aseModel','omegaMax');

%% LOAD STUFF AND START FROM HERE
if(~exist('ssObjs','var'))
    load('TIME_compare.mat')
end

%% PLOT COMPARE EXPERIMENTAL RESPONSE TO STATE-SPACE REPLICATION

for idxSpeed = 1:6
    for idxInput = 1:4
        plotTime(ssObjs(idxInput,idxSpeed),expObjs(idxInput,idxSpeed))
    end
end

%% SAVE PLOTS

wantSave = input("SAVE PLOTS? Y/N: ",'s');
if(wantSave=='Y')
    saveName = input("TYPE THE PLOT NAMING SUFFIX FOR SAVING THIS MODEL: ",'s');

    % create directory
    dirPath = ['./responsePlots/',saveName,'/'];
    if(~isfolder(dirPath))
        mkdir(dirPath);
    end

    % save in directory
    for idxSpeed = 1:6
        for idxInput = 1:4
            figure((idxSpeed-1)*4+idxInput)
            print([dirPath,'TIMECOMPARE_',saveName,'_',char(expObjs(idxInput,idxSpeed).fullTitle),'.png'],'-dpng','-r300')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = ssResponse(expObj,A,B,C,D)
    % import experimental conditions and input
    t_ = expObj.t;
    u_ = expObj.u;

    % % remove duplicate time records in experimental data
    % [t,mask] = unique(t_);
    % u = u_(mask,:);

    % downsample experimental data
    nStepSample = 20;
    t = t_(1:nStepSample:end);
    u = u_(1:nStepSample:end,:);
    
    % control input function
    uFunc = @(tq) interp1(t,u,tq,'makima')';

    % replace time vector
    t = [t(1),t(end)]; % specify endpoints only
    
    % derivative function
    odeFunc = @(t,x) A*x+B*uFunc(t);
    
    % zero initial state
    x0 = zeros(size(A,1),1);
    
    % integrate
    [tOut,xOut] = ode45(odeFunc,t,x0);
    xOut = xOut';
    
    % control history
    uOut = uFunc(tOut');
    
    % output history
    yOut = C*xOut+D*uOut;

    % return object
    obj.t = tOut;
    obj.u = uOut';
    obj.y = yOut';
end

%%
function plotTime(ssObj,expObj)
% INPUTS:
    % objs with fields t,u,y

    objs = {ssObj,expObj};
    labels = ["model","experiment"];
    fig = figure;
    fig.Position = [0,0,1000,750];
    tl = tiledlayout(6,1);
    xlabel(tl,'Time (s)')
    title(tl,expObj.fullTitle,'interpreter','none')

    % I/O names
    inputNames = ["ail1 (deg)","ail2 (deg)","elev (deg)","gust (deg)"];
    outputNames = ["acc1 (m/s2)","acc2 (m/s2)","acc3 (m/s2)","microstrain","pitch (deg)","pitchDot (deg/s)"];

    for idxObj = 1:length(objs)
        obj = objs{idxObj};

        t = obj.t;
        u = obj.u;
        y = obj.y;

        % plot input of experiment
        if(idxObj==1)
            ax = nexttile(1);
            plot(t,u)
            legend(inputNames)
            ylabel(ax,'Input')
            grid on
            hold on
        end
        
        % plot output
        for idxOut = 1:5 % acc1, acc2, acc3, strain, pitch (excluding pitchDot)
            ax = nexttile(idxOut+1);
            plot(t,y(:,idxOut),'DisplayName',labels(idxObj))
            ylabel(outputNames(idxOut))
            legend()
            grid on
            hold on
        end
    end
end