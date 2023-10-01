% Title: Plot MARGE NSTRAN Nodes
% Author: Anthony Su
% Date: 2023-05-14

close all
clear all
clc

% import node locations and convert to inches
xyz = load("ASEInputData/GRID_ID_XYZ.mat").GRID_ID_XYZ;
xyz(:,2:4) = xyz(:,2:4)*100/2.54;

% isolate wing geometry
wing_xyz = xyz(1:end-37,:); % wing

% plot wing
figure
plot(wing_xyz(:,2),wing_xyz(:,3),'k.','MarkerSize',10)
hold on
axis equal
xline(0)
yline(0)

for idx = 1:length(wing_xyz)
    id = wing_xyz(idx,1);
    x_ = wing_xyz(idx,2);
    y_ = wing_xyz(idx,3);

    % label coords
    if(id<2000) % if fuselage
        dy = -0.2;
        label = num2str(x_,'%.2f'); % label fuselage/root nodes with x-pos
    elseif(floor(mod(id/100,10))==4) % if x4xx
        dy = +0.4;
        % label = num2str(x_,'%.2f'); % label wing/tail nodes with x-pos
        label = num2str(y_,'%.2f'); % label wing/tail nodes with y-pos
    else
        dy = +0.2;
        % label = num2str(x_,'%.2f'); % label wing/tail nodes with x-pos
        label = num2str(y_,'%.2f'); % label wing/tail nodes with y-pos
    end
    text(x_,y_-dy,label,'FontSize',6)
    % label node id#
    label = num2str(id);
    text(x_,y_+1.5*dy,label,'FontSize',6)
end

ax = gca;
set(ax,'position',[0 0 1 1])
set(ax,'XTick',[])
set(ax,'YTick',[])
set(ax,'XLim',ax.XLim+[-2,2])


% plot 3D in WT
% figure
% plot3(xyz(:,2),xyz(:,3),xyz(:,4),'b.','MarkerSize',10) % WT
% axis equal
% grid on