% Title: Plot State-Space and Experimental FRFs Together
% Author: Anthony Su
% Date: 2023-08-17

% close all
clear all

d0 = load('FRF_ASE_SS.mat').dataObjs(1);
d1 = load('windTunnel\FRF_q207_A1S5.mat').dataObj;
d2 = load('windTunnel\FRF_q207_A2S5.mat').dataObj;
d3 = load('windTunnel\FRF_q207_ES2.mat').dataObj;
d4 = load('windTunnel\FRF_DG4.mat').dataObj;
dataObjs = [d0,d1,d2,d3,d4];
isExperiment = [0,1,1,1,1];

fileTitles = ["State-Space","q207_A1S5","q207_A2S5","q207_ES2"];

nInputs = 4;
nOutputs = 6;

outputNames = ["acc1 (g)","acc2 (g)","acc3 (g)","microstrain","pitch (deg)","pitch dot (deg/s)"];
inputNames = ["ail1 (deg)","ail2 (deg)","elev (deg)","gust vane (deg)"];

fig = figure;
fig.Position(1) = 0;
fig.Position(2) = 0;
fig.Position(3) = 1500;
fig.Position(4) = 750;
tiledlayout(nOutputs,nInputs);
for idxIn = 1:nInputs
    for idxOut = 1:nOutputs
        ax = nexttile((idxOut-1)*nInputs+idxIn);
        for idxObj = 1:length(dataObjs)
            obj = dataObjs(idxObj);
            freq = obj.freq(:,idxIn,idxOut);
            mag1 = abs(obj.H1_FRF(:,idxIn,idxOut)); % use H1 FRF because the input is noise-free (since it is the original commanded value)
            mag2 = abs(obj.H2_FRF(:,idxIn,idxOut)); % H2
            magv = abs(obj.Hv_FRF(:,idxIn,idxOut)); % Hv

            % magnitude to dB
            % mag1 = mag2db(mag1);
            % mag2 = mag2db(mag2);
            % magv = mag2db(magv);

            if(isExperiment(idxObj))
                plot(freq,mag1,'r-')
                hold on
                % plot(freq,mag2,'b-')
                % plot(freq,magv,'m-')
                grid on
            else
                plot(freq,magv,'k-')
                hold on
            end

            if(idxIn==1)
                ylabel(outputNames(idxOut))
            end
            if(idxOut==1)
                title(inputNames(idxIn))
            end
            if(idxOut~=nOutputs)
                set(ax,'xTickLabel',{})
            end
        end
    end
end

nexttile(4)
% legend(fileTitles,'Interpreter','none')
legend({'state-space','experiment H1 FRF','','','','experiment H2 FRF','','','','experiment Hv FRF','','',''})

tl = gcf().Children(1);
title(tl,'u=19 m/s','FontSize',24)
ylabel(tl,'Magnitude','FontSize',18)
% ylabel(tl,'Magnitude (dB)','FontSize',18)
xlabel(tl,'Frequency (Hz)','FontSize',18)