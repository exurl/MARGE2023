% Title: Plot State-Space and Experimental FRFs Together
% Author: Anthony Su
% Date: 2023-08-17

% loads experimental FRFs and state-space model FRFs and plots comparison

% close all
clear all

% speed indices of interest: pick which speeds to plot
speedIdxsInterest = [1:6];

% plot acceleromters?
global plotAcc
plotAcc = true;

%% IMPORT EXPERIMENTAL DATA

% input validation
assert(all(speedIdxsInterest>=1) & all(speedIdxsInterest<=6) & all(round(speedIdxsInterest)==speedIdxsInterest))

% speeds
q = [60,100,164,207,281,343];

% q60
if(any(speedIdxsInterest==1))
expObjs(1,1) = load('windTunnel\wtData\FRF_q60_A1S5.mat').data;
expObjs(2,1) = load('windTunnel\wtData\FRF_q60_A2S5.mat').data;
expObjs(3,1) = load('windTunnel\wtData\FRF_q60_ES2.mat').data;
expObjs(4,1) = load('windTunnel\wtData\FRF_q60_GD4.mat').data;
end

% q100
if(any(speedIdxsInterest==2))
expObjs(1,2) = load('windTunnel\wtData\FRF_q100_A1S5.mat').data;
expObjs(2,2) = load('windTunnel\wtData\FRF_q100_A2S5.mat').data;
expObjs(3,2) = load('windTunnel\wtData\FRF_q100_ES2.mat').data;
expObjs(4,2) = load('windTunnel\wtData\FRF_q100_GD4.mat').data;
end

% q164
if(any(speedIdxsInterest==3))
expObjs(1,3) = load('windTunnel\wtData\FRF_q164_A1S5.mat').data;
expObjs(2,3) = load('windTunnel\wtData\FRF_q164_A2S5.mat').data;
expObjs(3,3) = load('windTunnel\wtData\FRF_q164_ES2.mat').data;
expObjs(4,3) = load('windTunnel\wtData\FRF_q164_GD4.mat').data;
end

% q207
if(any(speedIdxsInterest==4))
expObjs(1,4) = load('windTunnel\wtData\FRF_q207_A1S5.mat').data;
expObjs(2,4) = load('windTunnel\wtData\FRF_q207_A2S5.mat').data;
expObjs(3,4) = load('windTunnel\wtData\FRF_q207_ES2.mat').data;
expObjs(4,4) = load('windTunnel\wtData\FRF_q207_GD4.mat').data;
end

% q281
if(any(speedIdxsInterest==5))
expObjs(1,5) = load('windTunnel\wtData\FRF_q281_A1S5.mat').data;
expObjs(2,5) = load('windTunnel\wtData\FRF_q281_A2S5.mat').data;
expObjs(3,5) = load('windTunnel\wtData\FRF_q281_ES2.mat').data;
expObjs(4,5) = load('windTunnel\wtData\FRF_q281_GD4.mat').data;
end

% q343
if(any(speedIdxsInterest==6))
expObjs(1,6) = load('windTunnel\wtData\FRF_q343_A1S3p5.mat').data;
expObjs(2,6) = load('windTunnel\wtData\FRF_q343_A2S3p5.mat').data;
expObjs(3,6) = load('windTunnel\wtData\FRF_q343_ES1.mat').data;
expObjs(4,6) = load('windTunnel\wtData\FRF_q343_GD4.mat').data;
end

% take only relevant FRFs from each experiment
for idxSpeed = speedIdxsInterest
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
for idxSpeed = speedIdxsInterest
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

% ask to save plots
wantSave = input("SAVE PLOTS? Y/N: ",'s');
if(wantSave=='Y')
    saveName = input("TYPE THE PLOT NAMING SUFFIX FOR SAVING THIS MODEL: ",'s');

    % create directory
    dirPath = ['./responsePlots/',saveName,'/'];
    if(~isfolder(dirPath))
        mkdir(dirPath);
    end

    % set plot visibility
    visibility = 'off';
else
    visibility = 'on';
end

% plot for each speed
for idxSpeed = speedIdxsInterest
    % plot
    plotFreq(ssObjs(:,idxSpeed),expObjs(:,idxSpeed),q(idxSpeed),visibility)

    % save in directory
    if(wantSave=='Y')
        savename = [dirPath,'FRFCOMPARE_',saveName,'_q',char(num2str(q(idxSpeed))),'.png'];
        print(savename,'-dpng','-r300')
        disp(['saved ',savename])
    end
end

% save plots
% wantSave = input("SAVE PLOTS? Y/N: ",'s');
% if(wantSave=='Y')
%     saveName = input("TYPE THE PLOT NAMING SUFFIX FOR SAVING THIS MODEL: ",'s');
% 
%     % get only relevant plots
%     rootObj = groot;
%     figs = [rootObj.Children];
%     for i = 1:length(figs)
%         try
%             figMask(i) = isequal(figs(i).Children(1).Children(1).String,{'experiment Hv FRF','model'});
%         catch
%             figMask(i) = false;
%         end
%     end
%     figNumbers = [figs(figMask).Number];
%     if(length(figNumbers)>6)
%         error('too many relevant plots open; not sure which ones to save!')
%     end
%     figNumbers = sort(figNumbers);
% 
%     % create directory
%     dirPath = ['./responsePlots/',saveName,'/'];
%     if(~isfolder(dirPath))
%         mkdir(dirPath);
%     end
% 
%     % save in directory
%     for idxSpeed = 1:length(speedIdxsInterest)
%         figure(figNumbers(idxSpeed))
%         savename = [dirPath,'FRFCOMPARE_',saveName,'_q',char(num2str(q(idxSpeed))),'.png'];
%         print(savename,'-dpng','-r300')
%         disp(['saved ',savename])
%     end
% end


%%
%% FUCTIONS
%%

function plotFreq(ssObjs,expObjs,q,visibility)
% INPUTS: ssObjs and expObjs have fields freq, H1_FRF, H2_FRF, Hv_FRF
% they have dims (nInputs) and internally the fields have dims (:,nOutputs)

    global plotAcc

    % define utility variables
    nInputs = 4;
    nOutputs = 6;
    outputNames = ["acc1 (g)","acc2 (g)","acc3 (g)","microstrain","pitch (deg)","pitch dot (deg/s)"];
    inputNames = ["ail1 (deg)","ail2 (deg)","elev (deg)","gust vane (deg)"];
    
    % initialize plot
    fig = figure('visible',visibility);
    fig.Position(1) = 0;
    fig.Position(2) = 0;
    fig.Position(3) = 1920;
    fig.Position(4) = 1080;

    % plot accels?
    if(plotAcc==false)
        idxOut0 = 4;
        accShift = -3;
        tl = tiledlayout(nOutputs-3,nInputs,'Padding', 'none', 'TileSpacing', 'compact');
    else
        idxOut0 = 1;
        accShift = 0;
        tl = tiledlayout(nOutputs,nInputs,'Padding', 'none', 'TileSpacing', 'compact');
    end
    
    % plot data
    for idxIn = 1:nInputs
        expObj = expObjs(idxIn);
        ssObj = ssObjs(idxIn);
        for idxOut = idxOut0:nOutputs
            % initialize axes
            ax = nexttile((idxOut+accShift-1)*nInputs+idxIn);
            
            % plot experiment
            % plot(expObj.freq(:,idxOut),abs(expObj.Hv_FRF(:,idxOut)),'-k','DisplayName','experiment')
            b_ = [abs(expObj.Hv_FRF(:,idxOut)-expObj.H1_FRF(:,idxOut)),abs(expObj.Hv_FRF(:,idxOut)-expObj.H2_FRF(:,idxOut))];
            [lineObj_,patchObj_] = boundedline(expObj.freq(:,idxOut),abs(expObj.Hv_FRF(:,idxOut)),b_,'-k');
            set(lineObj_,'DisplayName','experiment Hv FRF')
            set(patchObj_,'HandleVisibility','off')
            hold on
            % plot(expObj.freq(:,idxOut),abs(expObj.H1_FRF(:,idxOut)),'--k','DisplayName','experiment H1')
            % plot(expObj.freq(:,idxOut),abs(expObj.H2_FRF(:,idxOut)),':k','DisplayName','experiment H2')

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
                % title(ax,expObj.title,inputNames(idxIn))
                title(ax,'',inputNames(idxIn))
            end

            % remove excess ticks, labels
            if(idxOut<nOutputs)
                ax.XTickLabels = [];
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

    % grey out acc3 and pitchDot outputs which are unreliable
    % for idxIn = 1:nInputs
    %     for idxOut = [3,6]
    %         ax = nexttile((idxOut-1)*nInputs+idxIn);
    %         ax.Color = [1,0.9,0.9];
    %     end
    % end

end