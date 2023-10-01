% Title: ns vs.Strain Study
% Author: Anthony Su
% Date: 2023-08-24

clear all

% mode numbers of interest
nsVec = 2:15;

% import models
for idxFile = 1:length(nsVec)
    ns = nsVec(idxFile);
    aseModels(idxFile) = load(strcat('ASE_SS_ns',num2str(ns),'.mat'));
end

% initialize I/O names
inputNames = ["aileron 1","aileron 2","elevator","gust vane"];
outputNames = ["wing LE acc","wing TE acc","tail acc","root strain","pitch","pitch first derivative"];

omega = linspace(0.5,4,1e4); % Hz
nc = 4;
magnitudeResponses = zeros(length(omega),nc,length(aseModels));

for idxFile = 1:length(nsVec)
    % bring model data into workspace
    aseModel = aseModels(idxFile);
    A = aseModel.A;
    Bc = aseModel.Bc;
    C = aseModel.C;
    Dc = aseModel.Dc;
    omegaMax = aseModel.omegaMax;
    rho = aseModel.rho;
    u = aseModel.u';
    q = 0.5*rho*u.*u;
    nSpeeds = length(u);
    nStates = size(A,1);
    nInputs = size(Bc,2);
    nOutputs = size(C,1);

    % take only strain output
    C = C(4,:,:);
    Dc = Dc(4,:,:);

    % frequency response
    sys = ss(A,Bc,C,Dc);
    for idxSpeed = 1:nSpeeds
        sys_ = sys(:,:,idxSpeed);
        [mag(:,:,:,idxSpeed),phase(:,:,:,idxSpeed),] = bode(sys_,omega*2*pi);
            % indices of mag and phase are (idxOutput,idxInput,:,idxSpeed)
    end
    
    % take only the 19 m/s data
    mag = mag(:,:,:,1);
    phase = phase(:,:,:,1);

    % reshape
    mag = squeeze(mag)';     % new shape: (:,idxInput)
    phase = squeeze(phase)'; % new shape: (:,idxInput)

    % store data
    magnitudeResponses(:,:,idxFile) = mag;
    clear mag phase
end

% convert from absolute magnitude to dB
% magnitudeResponses = mag2db(magnitudeResponses);

%% PLOT RESPONSES
styles = ["r-","m-","b-","k:","k:","k:","k:","k:","k:","k:","k:","k:","k:","k:"];

fig = figure;
fig.Position = [0,0,1600,1000];

% meta-tiledlayout
TL = tiledlayout(1,3);
ylabel(TL,'Magnitude','FontSize',14)
title(TL,'Strain Gauge Response to Actuators with Various Fidelity Models','FontSize',24)

tl = tiledlayout(TL,4,1);
tl.Layout.Tile = 1;
title(tl,'Response','FontSize',18)
xlabel(tl,'Frequency (Hz)','FontSize',14)

% plot input
for idxInput = 1:nc
    ax = nexttile(tl,idxInput);
    hold on
    for idxFile = 1:length(nsVec)
        ns = nsVec(idxFile);
        plot(omega,magnitudeResponses(:,idxInput,idxFile),styles(idxFile),'DisplayName',strcat('ns=',num2str(ns)))
    end
    title(ax,inputNames(idxInput))
    grid on

    % remove x-tick labels except for bottom plot
    if(idxInput~=nc)
        set(ax,'XTickLabel',{})
    end
end
% leg = legend;
leg.Layout.Tile = 'east';

%% PLOT DIFFERENCE BETWEEN RESPONSES AND NS=15
% fig = figure;
% fig.Position(1:2) = [0,0];
% fig.Position(3:4) = fig.Position(3:4)*1.2;
tl = tiledlayout(TL,4,1);
tl.Layout.Tile = 2;
title(tl,'Error (w/ ns=15 as truth)','FontSize',18)
xlabel(tl,'Frequency (Hz)','FontSize',14)

% plot input
for idxInput = 1:nc
    ax = nexttile(tl,idxInput);
    hold on
    for idxFile = 1:length(nsVec)-1
        ns = nsVec(idxFile);
        magDiff = magnitudeResponses(:,idxInput,length(nsVec))-magnitudeResponses(:,idxInput,idxFile);
        plot(omega,magDiff,styles(idxFile),'DisplayName',strcat('ns=',num2str(ns)))
    end
    title(ax,inputNames(idxInput))
    grid on

    % remove x-tick labels except for bottom plot
    if(idxInput~=nc)
        set(ax,'XTickLabel',{})
    end
end
% leg = legend;
leg.Layout.Tile = 'east';

%% PLOT RATIO BETWEEN RESPONSES AND NS=15
% fig = figure;
% fig.Position(1:2) = [0,0];
% fig.Position(3:4) = fig.Position(3:4)*1.2;
tl = tiledlayout(TL,4,1);
tl.Layout.Tile = 3;
title(tl,'Error Fraction (w/ ns=15 as truth)','FontSize',18)
xlabel(tl,'Frequency (Hz)','FontSize',14)

% plot input
for idxInput = 1:nc
    ax = nexttile(tl,idxInput);
    hold on
    for idxFile = 1:length(nsVec)-1
        ns = nsVec(idxFile);
        magDiff = magnitudeResponses(:,idxInput,length(nsVec))-magnitudeResponses(:,idxInput,idxFile);
        ratio = magDiff./magnitudeResponses(:,idxInput,length(nsVec));
        plot(omega,ratio,styles(idxFile),'DisplayName',strcat('ns=',num2str(ns)))
    end
    title(ax,inputNames(idxInput))
    grid on

    % remove x-tick labels except for bottom plot
    if(idxInput~=nc)
        set(ax,'XTickLabel',{})
    end
end
leg = legend;
leg.Layout.Tile = 'east';