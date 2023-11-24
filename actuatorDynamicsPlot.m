% Title: Actuator Dynamics Plots
% Author: Anthony Su
% Date: 2023-11-23

close all
clear all
clc

%% gust vanes
fig = figure;
pade(0.34,2);

% adjust title
fig.Children(1).Title.String = 'Phase response';
fig.Children(3).Title.String = 'Step response';

% add grid
grid(fig.Children(1),'on')
grid(fig.Children(3),'on')

% adjust font size
fig.Children(1).FontSize = 9.9;
fig.Children(3).FontSize = 9.9;

% frequency axis units from rad to Hz
fig.Children(1).Children(1).XData = rad2deg(fig.Children(1).Children(1).XData);
fig.Children(1).Children(2).XData = rad2deg(fig.Children(1).Children(2).XData);
xlim([fig.Children(1).Children(1).XData(1),fig.Children(1).Children(1).XData(end)])
fig.Children(1).XLabel.String = 'Frequency (Hz)';

% export plot
wantSave = input('SAVE GUST VANE DYNAMICS PLOT? Y/N: ','s');
if(wantSave == 'Y')
    print('gustInputDynamics','-dpng','-r300')
end

%% servos
G = tf([0,0,1461],[1,62.2,1461]);


fig = figure;
bodeOpt = bodeoptions;

% adjust frequnecy units
bodeOpt.FreqUnits = 'Hz';

% adjust font size
bodeOpt.XLabel.FontSize = 9.9;
bodeOpt.YLabel.FontSize = 9.9;

% plot
bodeplot(G,bodeOpt)

% adjust font size again
fig.Children(2).FontSize = 9.9;
fig.Children(3).FontSize = 9.9;

% add grid
grid(fig.Children(2),'on')
grid(fig.Children(3),'on')

% remove title
title('')

% export plot
wantSave = input('SAVE SERVO DYNAMICS PLOT? Y/N: ','s');
if(wantSave == 'Y')
    print('servoInputDynamics','-dpng','-r300')
end