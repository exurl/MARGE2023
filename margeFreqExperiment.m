% Title: Plot State-Space and Experimental FRFs Together
% Author: Anthony Su
% Date: 2023-08-17

close all
clear all

%% IMPORT EXPERIMENTAL DATA

% speeds
q = [60,100,164,207,281,343];

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

% take only relevant FRFs from each experiment
for idxSpeed = 1:6
for idxInput = 1:4
    obj.freq = squeeze(expObjs(idxInput,idxSpeed).freq(:,idxInput,:));
    obj.H1_FRF = squeeze(expObjs(idxInput,idxSpeed).H1_FRF(:,idxInput,:));
    obj.H2_FRF = squeeze(expObjs(idxInput,idxSpeed).H2_FRF(:,idxInput,:));
    obj.Hv_FRF = squeeze(expObjs(idxInput,idxSpeed).Hv_FRF(:,idxInput,:));
    obj.title = expObjs(idxInput,idxSpeed).title;
    newExpObjs(idxInput,idxSpeed) = obj;
end
end
expObjs = newExpObjs;
clear newExpObjs

%% IMPORT STATE-SPACE MODEL DATA

% import
ssSuperObjs = load('FRF_ASE_SS.mat').dataObjs; % (1)x(6) array corresponding to speeds

% break each obj into 4 inputs
for idxSpeed = 1:6
    for idxInput = 1:4
        obj.freq = squeeze(ssSuperObjs(idxSpeed).freq(:,idxInput,:));
        obj.H1_FRF = squeeze(ssSuperObjs(idxSpeed).H1_FRF(:,idxInput,:));
        obj.H2_FRF = squeeze(ssSuperObjs(idxSpeed).H2_FRF(:,idxInput,:));
        obj.Hv_FRF = squeeze(ssSuperObjs(idxSpeed).Hv_FRF(:,idxInput,:));
        obj.title = ssSuperObjs(idxSpeed).title;
        ssObjs(idxInput,idxSpeed) = obj;
    end
end


%% PLOT

% plot for each speed
for idxSpeed = 1:6
% for idxSpeed = 6 % plot q343 only
    plotFreq(ssObjs(:,idxSpeed),expObjs(:,idxSpeed),q(idxSpeed))
end

% save plots
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
        figure(idxSpeed)
        print([dirPath,'FRFCOMPARE_',saveName,'_q',char(num2str(q(idxSpeed))),'.png'],'-dpng','-r300')
    end
end


%%
%% FUCTIONS
%%

function plotFreq(ssObjs,expObjs,q)
% INPUTS: ssObjs and expObjs have fields freq, H1_FRF, H2_FRF, Hv_FRF
% they have dims (nInputs) and internally the fields have dims (:,nOutputs)

    % define utility variables
    nInputs = 4;
    nOutputs = 6;
    outputNames = ["acc1 (g)","acc2 (g)","acc3 (g)","microstrain","pitch (deg)","pitch dot (deg/s)"];
    inputNames = ["ail1 (deg)","ail2 (deg)","elev (deg)","gust vane (deg)"];
    
    % initialize plot
    fig = figure;
    fig.Position(1) = 0;
    fig.Position(2) = 0;
    fig.Position(3) = 1920;
    fig.Position(4) = 1080;
    tl = tiledlayout(nOutputs,nInputs);

    % plot data
    for idxIn = 1:nInputs
        expObj = expObjs(idxIn);
        ssObj = ssObjs(idxIn);
        for idxOut = 1:nOutputs
            % initialize axes
            ax = nexttile((idxOut-1)*nInputs+idxIn);
            
            % plot experiment
            plot(expObj.freq(:,idxOut),abs(expObj.Hv_FRF(:,idxOut)),'-k','DisplayName','experiment')
            hold on

            % plot state-space model
            plot(ssObj.freq(:,idxOut),abs(ssObj.Hv_FRF(:,idxOut)),'-r','DisplayName','model')

            % formatting
            grid on

            % y-axis output label
            if(idxIn==1)
                ylabel(outputNames(idxOut))
            end

            % x-axis input label
            if (idxOut==1)
                title(ax,expObj.title,inputNames(idxIn))
            end
        end
    end

    % legend
    lg  = legend(); 
    lg.Layout.Tile = 'North'; % place to the right of the tiledlayout
    
    % formatting
    title(tl,['q=',char(num2str(q))],'FontSize',24)
    ylabel(tl,'Magnitude','FontSize',18)
    xlabel(tl,'Frequency (Hz)','FontSize',18)
    warning('off',char([])) % remove warning about linkaxes in next line
    linkaxes(tl.Children,'x') % unify x-axes
end