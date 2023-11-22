% Title: Plot MARGE NASTRAN Nodes of Accelerometers in GVT
% Author: Anthony Su
% Date: 2023-05-14 (modified 2023-07/18)

close all
clear all
clc

% import node locations
xyz = load("ASEInputData/GRID_ID_XYZ.mat").GRID_ID_XYZ;

% isolate wing geometry
wing_xyz = xyz(1:end-37,:); % wing

% append fake nodes
spar_width = 3/100*2.54; % 3 inch --> meters
spar_LE = wing_xyz(wing_xyz(:,1)==2202,2);
spar_TE = spar_LE+spar_width;

node2314 = wing_xyz(wing_xyz(:,1)==2214,:);
node2314(1) = 2314;
node2314(2) = spar_TE;
wing_xyz = [wing_xyz;node2314];

node2224 = wing_xyz(wing_xyz(:,1)==2024,:);
node2224(1) = 2224;
node2224(2) = spar_LE;
wing_xyz = [wing_xyz;node2224];

node2324 = wing_xyz(wing_xyz(:,1)==2024,:);
node2324(1) = 2324;
node2324(2) = spar_TE;
wing_xyz = [wing_xyz;node2324];

node3208 = wing_xyz(wing_xyz(:,1)==3208,:);
node3308 = node3208;
node3308(1) = 3308;
node3308(2) = node3208(2)+3;
wing_xyz = [wing_xyz;node3308];

node2202 = wing_xyz(wing_xyz(:,1)==2202,:);
nodeRoot = node2202;
nodeRoot(1) = 4001;
nodeRoot(2) = 0;
wing_xyz = [wing_xyz;nodeRoot];

node2223 = wing_xyz(wing_xyz(:,1)==2223,:);
nodeTip = node2223;
nodeTip(1) = 4002;
nodeTip(2) = 0;
wing_xyz = [wing_xyz;nodeTip];

nodeMassTip = node2224;
nodeMassTip(1) = 4003;
nodeMassTip(2) = nodeMassTip(2)-0.5/100*2.54;
wing_xyz = [wing_xyz;nodeMassTip];

nodeMassTail = node2324;
nodeMassTail(1) = 4004;
nodeMassTail(2) = nodeMassTail(2)+0.5/100*2.54;
wing_xyz = [wing_xyz;nodeMassTail];

% plot wing
fig = figure;
fig.Position = [10,10,800,600];
border_ids = [1001,2001,2002,4001,4002,2023,2024,4003,4004,2024,2023,2823,2802,2002,2001,3001,3002,3202,3208,3808,3802,3002];
x_ = [];
y_ = [];
for idx = 1:length(border_ids)
    node = wing_xyz(wing_xyz(:,1)==border_ids(idx),:);
    x_ = [x_,node(2)];
    y_ = [y_,node(3)];
end
plot(x_,y_,'k-')
hold on
ax = gca;
axis equal
xlim([-0.1,0.9])
ylim([-0.1,0.6])

% isolate tested nodes
GVT1_ids = sort([2014,2020,2024,1017,1034]);
GVT2_ids = sort([2012,2214,2314,2020,2224,2324,1013,1028,3208,3308]);

% plot tested nodes
GVT1_nodes = [];
for idx = 1:length(GVT1_ids)
    node = wing_xyz(wing_xyz(:,1)==GVT1_ids(idx),:);
    GVT1_nodes = [GVT1_nodes;node];
end
GVT2_nodes = [];
for idx = 1:length(GVT2_ids)
    node = wing_xyz(wing_xyz(:,1)==GVT2_ids(idx),:);
    GVT2_nodes = [GVT2_nodes;node];
end

plot(GVT1_nodes(:,2),GVT1_nodes(:,3),'bo','MarkerSize',6,'DisplayName','dataset 1')
plot(GVT2_nodes(:,2),GVT2_nodes(:,3),'r*','MarkerSize',6,'DisplayName','dataset 2')

% plot tested node names
rotation = 30;
for idx = 1:length(GVT1_ids)
    name_ = ['  ',char(string(GVT1_nodes(idx,1)))];
    x_ = GVT1_nodes(idx,2);
    y_ = GVT1_nodes(idx,3);
    txt_ = text(x_,y_,name_);
    txt_.Rotation = rotation;
end
for idx = 1:length(GVT2_ids)
    name_ = ['  ',char(string(GVT2_nodes(idx,1)))];
    x_ = GVT2_nodes(idx,2);
    y_ = GVT2_nodes(idx,3);
    txt_ = text(x_,y_,name_);
    txt_.Rotation = rotation;
end

ax = gca;
% set(ax,'position',[0 0 1 1])
% set(ax,'XTick',[])
% set(ax,'YTick',[])

% flip y-axis
set(ax,'ydir','reverse')

legend([ax.Children(end-1),ax.Children(end-2)],'Location','southeast')
xlabel('meters')
ylabel('meters')

% export with tight margins
% exportgraphics(ax,'accelLocPlot.png','Resolution',300)