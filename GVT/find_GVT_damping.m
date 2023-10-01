% Title: Find Damping
% Author: Anthony Su
% Date: 2022-07-11

% to use this code, adjust the following variables:
    % baseDir
    % selectLocs
    % windows

% close all
clear all
clc

%% PROGRAM PARAMETERS

% directory to search for FRF data files
% baseDir = 'GVT1_reformatted/'
baseDir = 'GVT2/'
% baseDir = './'

% select files of location of interest
% selectLocs = ["1017z","1034z","3008z","3408z"] % GVT1 fuse bend nodes
% selectLocs = ["2012z","2214z","2314z","2020z","2224z","2324z"] % GVT2 wing
% selectLocs = ["2014z","2020z","2024z"]; % GVT1 wing bending
% selectLocs = ["1034y","2024x"] % GVT1 in-plane
selectLocs = ["2214x","2224x"] % GVT2 in-plane

% windows containing modes of interest
% windows = [
%     19,20.5; % fuse bend 1
%     71,77] % fuse bend 2
% windows = [8,12; % wing bend 2
%     13,21; % wing twist 1
%     31.5,34; wing bend 3
%     40,57; % wing twist 2
%     56,65] % wing bend 4
% windows = [1.3,1.55]; % wing bend 1
% windows = [18,22] % fuse in-plane bend 1
windows = [60,80]; % wing in-plane bending 1 <-- shotgun blast wide guess

% names of modes of interest corresponding to windows
% modeNames = ["1st fuselage bending","2nd fuselage bending"]
% modeNames = ["wing bending 2","wing torsion 1","wing bending 3","wing torsion 2","wing bending 4"]
% modeNames = ["wing bending 1"]
% modeNames = ["fuselage in-plane bending 1"]
modeNames = ["wing in-plane bending 1"]

%% IMPORT DATA

% find files
files = dir(baseDir);
fileNames = strings(0);
for idxFile = 1:length(files)
    f_ = files(idxFile).name;
    if(contains(f_,'FRF_') && contains(f_,'.mat'))
    if(contains(f_,selectLocs))
        fileNames(end+1) = string(f_);
    end
    end
end
N_files = length(fileNames);
filePaths = strings(size(fileNames));
for idxFile = 1:N_files
    filePaths(idxFile) = string([baseDir,char(fileNames(idxFile))]);
end

% check files exist
if(any(~isfile(filePaths)))   
    error('file does not exist!')
end

% import data
for idxFile = 1:N_files
    % import data object
    f = load(filePaths(idxFile));
    if(idxFile==1)
        % initialize array
        FRFs(1) = f;
    else
        % append to array
        FRFs(end+1) = f;
    end
end

%% COMPUTE DAMPING RATIOS

% compute damping, natural frequency, half-power frequencies
[zeta,fMax,f1,f2] = getDamping(FRFs,windows);

% reshape data for 2-D format
FRFTypeNames = ["H1","H2","Hv"];
for idxMode = 1:length(modeNames)
    disp('=============================================')
    disp(['MODE: ',char(modeNames(idxMode))])
    disp('=============================================')
    
    for idxFRFType = [3,1,2] % show Hv first
        disp(['-- ',char(FRFTypeNames(idxFRFType)),' ---------------------------------------'])
        disp('  ___fMax     _____f1     _____f2      __zeta')
        matOut = [reshape(fMax(idxMode,1,:),[],1),...
            reshape(f1(idxMode,idxFRFType,:),[],1),...
            reshape(f2(idxMode,idxFRFType,:),[],1),...
            reshape(zeta(idxMode,idxFRFType,:),[],1)...
            ];
        printmatrix(matOut); % prints in format copy-pastable to excel
    end
end

% also print file origin of results
disp('=============================================')
disp('FILES OF ORIGIN')
disp('=============================================')
for f_ = filePaths
    disp(f_)
end

%% FIND WING-Z DAMPING RATIOS (OUTDATED)

% % windows to search, each containing one peak only
% windows = [8,12;
%     13,21;
%     31.5,34;
%     40,57;
%     56,65];
% mode_names = ["wing bend 2","wing twist 1","wing bend 3","wing twist 2","wing bend 4"];
% 
% % find damping
% [zeta,f_max,f_1,f_2] = getDamping(FRFs{1},windows); % data from wing z 1
% [z_,f_,f1_,f2_] = getDamping(FRFs{2},windows); % data from wing z 2
% zeta = [zeta,z_];
% f_max = [f_max,f_];
% f_1 = [f_1,f1_];
% f_2 = [f_2,f2_];
% [z_,f_,f1_,f2_] = getDamping(FRFs{2},windows); % data from twist
% zeta = [zeta,z_];
% f_max = [f_max,f_];
% f_1 = [f_1,f1_];
% f_2 = [f_2,f2_];
% 
% % plot damping
% plotDamping(zeta,windows,mode_names)

%% FIND FUSE-Z DAMPING RATIOS (OUTDATED)

% % windows to search, each containing one peak only
% windows = [1.5,4;
%     17,22;
%     64,80;
%     ];
% mode_names = ["fuse bend 1","fuse bend 2","fuse bend 3"];
% 
% % find damping
% [zeta,f_max,f_1,f_2] = getDamping(FRFs{3},windows); % data from fuse z
% 
% % plot damping
% plotDamping(zeta,windows,mode_names)

%% FIND IN-PLANE DAMPING RATIOS (OUTDATED)

% % windows to search, each containing one peak only
% windows = [16.5,20];
% mode_names = ["fuse in-plane 1"];
% 
% % find damping
% [zeta,f_max,f_1,f_2] = getDamping(FRFs{4},windows); % data from minor
% 
% % plot damping
% plotDamping(zeta,windows,mode_names)

%% FIND WINGZ-LONG DAMPING RATIOS (OUTDATED)

% % windows to search, each containing one peak only
% windows = [1.3,1.55];
% mode_names = ["wing bend 1"];
% 
% % find damping
% [zeta,f_max,f_1,f_2] = getDamping(FRFs{6},windows); % data from minor
% 
% % plot damping
% plotDamping(zeta,windows,mode_names,0.002)

%% FUCTIONS

function [zeta,fMax,f1,f2] = getDamping(FRFs,windows)
% INPUTS:
    % FRFs:    array of FRF objects
    % windows: N x 2 matrix, each row is a search interval
% OUTPUTS:
    % zeta:    damping of Hv, H1, and H2 FRFs for each mode
    % fMax:    natural frequency of Hv, H1, and H2 FRFs for each mode
    % f1:      lower half-power freq of Hv, H1, and H2 FRFs for each mode
    % f2:      upper half-power freq of Hv, H1, and H2 FRFs for each mode

    % for each wing-z FRF, for each peak, find damping
    zeta = zeros(size(windows,1),3,length(FRFs));
    fMax = zeros(size(windows,1),3,length(FRFs));
    f1 = zeros(size(windows,1),3,length(FRFs));
    f2 = zeros(size(windows,1),3,length(FRFs));
    
    % for each FRF
    for idxFRF = 1:length(FRFs)
        f = FRFs(idxFRF).freq;

        % for each type of FRF
        for idxFRFType = 1:3
            switch idxFRFType
                case 1
                    FRF = FRFs(idxFRF).mag_H1;
                case 2
                    FRF = FRFs(idxFRF).mag_H2;
                case 3
                    % H_v: "robust TF estimate", Rocklin et. al. 1985
                    FRF = 0.5*(FRFs(idxFRF).mag_H1+FRFs(idxFRF).mag_H2);
            end

            % for each search window
            for idxWindow = 1:size(windows,1)
                % constrain data to search window
                mask = f >= windows(idxWindow,1) & f <= windows(idxWindow,2);
            
                % find damping, natural freq, half-power freqs
                try
                    [z_,f_,f1_,f2_] = halfPowerDamping(f(mask),FRF(mask));
                    zeta(idxWindow,idxFRFType,idxFRF) = z_;
                    fMax(idxWindow,idxFRFType,idxFRF) = f_;
                    f1(idxWindow,idxFRFType,idxFRF) = f1_;
                    f2(idxWindow,idxFRFType,idxFRF) = f2_;
                catch
                    % not found in window --> save as NaN
                    zeta(idxWindow,idxFRFType,idxFRF) = NaN;
                    fMax(idxWindow,idxFRFType,idxFRF) = NaN;
                    f1(idxWindow,idxFRFType,idxFRF) = NaN;
                    f2(idxWindow,idxFRFType,idxFRF) = NaN;
                end
            end
        end
    end
end



function plotDamping(zeta,windows,mode_names,bin_width)
% creates histogram of damping values computed with median highlighted
    
    % default bin width if unspecified
    if(~exist('bin_width','var'))
        disp('bin_width defaulted to 0.005')
        bin_width = 0.005;
    end

    % create+format figure+tiledlayout
    f = figure;
    xlab = matlab.graphics.layout.Text;
    set(xlab,'String','\zeta');
    tl = tiledlayout('vertical','TileSpacing','tight','Padding','tight','XLabel',xlab);
    for idx_window = 1:size(windows,1)
        if(idx_window==1)
            axs = nexttile(idx_window);
        else
            ax = nexttile(idx_window);
            axs = [axs,ax];
        end
    
        % plot
        histogram(zeta(idx_window,:),'BinWidth',bin_width)
        hold on
        xline(median(zeta(idx_window,:),'omitmissing'),'-r','LineWidth',2)
    
        % label+format
        title(mode_names(idx_window))
        if(idx_window==1)
            legend({'','median'})
        end
        grid on
    end
    linkaxes(axs,'xy')
end