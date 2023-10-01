% Title: Import Wind Tunnel Data
% Author: Anthony Su
% Date: 2023-08-11

% this script is for testing purposes. The goal is to reuse part of the
% code in this script in an upcoming more integrated wind tunnel data
% postprocessing pipeline

% close all
clear all
clc

%% SELECT FILES

% files = ["2023 summer/20230719/8ms CG2/output.mat",...
    % "2023 summer/20230719/8ms CG4/output.mat",...
    % "2023 summer/20230719/8ms ES2/output.mat"];

files = ["2023 summer/20230816/q207_A1S5_R1/output.mat",
    "2023 summer/20230816/q207_A1S5_R2/output.mat",
    "2023 summer/20230816/q207_A1S5_R3/output.mat",
    "2023 summer/20230816/q207_A2S5_R1/output.mat",
    "2023 summer/20230816/q207_A2S5_R2/output.mat",
    "2023 summer/20230816/q207_A2S5_R3/output.mat",
    "2023 summer/20230816/q207_ES2_R1/output.mat",
    "2023 summer/20230816/q207_ES2_R2/output.mat",
    "2023 summer/20230816/q207_ES2_R3/output.mat"];
    % q=207 Pa --> ~u=19 m/s

% check files exist
if(any(~isfile(files)))   
    error('file does not exist!')
end

dataObjs = [];
for idxFile = 1:length(files)
    %% IMPORT
    
    % load file
    path = files(idxFile);
    f = load(path);
    d = f.data;

    % data name
    pathParts = split(files(idxFile),'/');
    mask = find(contains(pathParts,'_'));
    dataName = pathParts(mask);

    
    % find locations of relevant data in object
    elemNames = string(getElementNames(d));
    idxCmd = find(elemNames=="cmd");
    idxGust = find(elemNames=="gust cmd");
    idxPlant = find(elemNames=="Plant:1");
    idxAcc1 = find(elemNames=="accel1_filtered");
    idxAcc2 = find(elemNames=="accel2_filtered");
    idxAcc3 = find(elemNames=="accel3_filtered");
    idxStrain = find(elemNames=="microstrain");
    
    % input data
    u_elev = d{idxCmd}.Values.elevator.Data;
    u_ail1 = d{idxCmd}.Values.aileron_1.Data;
    u_ail2 = d{idxCmd}.Values.aileron_2.Data;
    u_gust = d{idxGust}.Values.Data;
    u = [u_ail1,u_ail2,u_elev,u_gust];
    
    % output data 
    theta = d{idxPlant}.Values.theta.Data;
    acc1 = d{idxAcc1}.Values.Data;
    acc2 = d{idxAcc2}.Values.Data;
    acc3 = d{idxAcc3}.Values.Data;
    strain = d{idxStrain}.Values.Data;
    y = [acc1,acc2,acc3,strain,theta];
    
    % time vector
    t = d{idxStrain}.Values.Time;
    rate = (length(t)-1)/(t(end)-t(1));

    %% COMPUTE DERIVATIVES

    % first and second derivative of all inputs
    u = [u,[diff(u,1,1);zeros(1,4)],[zeros(1,4);diff(u,2,1);zeros(1,4)]];
        % ^zeros appended/prepended to derivatives to maintain series length
    
    % first derivative of pitch angle output
    y = [y,[diff(theta)*rate;0]];


    %% SHIFT AND TRUCATE
    
    % % truncate before input starts, 2s after input ends
    inputMask = find(sum(u~=0,2)); % indices with nonzero input
    idx1 = inputMask(1)-1; % start 1 datapoint before input
    idx2 = inputMask(end)+2*rate; % end 2 s after input

    % truncate before input starts, after most active signal drops to 20% of max magnitude
    % inputMask = find(sum(u~=0,2)); % indices with nonzero input
    % idx1 = inputMask(1)-1; % start 1 datapoint before input
    % [~,idxActive] = max(var(y./max(y)));
    % yActive = y(:,idxActive);
    % activeIdxs = find(abs(yActive)>0.2*max(abs(yActive)));
    % idx2 = activeIdxs(end);

    if(idx2>length(y))
        idx2 = length(y);
    end
    if(idx1==0)
        idx1 = 1;
    end
    inputRange = idx1:idx2;
    t = t(inputRange);
    u = u(inputRange,:);
    y = y(inputRange,:);

    % set new start time to zero
    t = t-t(1);

    % zero-mean everything
    u = u-mean(u);
    y = y-mean(y);
    
    % store as object
    obj.path = path;
    obj.title = dataName;
    obj.u = u;
    obj.y = y;
    obj.t = t;
    obj.rate = rate;
    dataObjs = [dataObjs,obj];
    clear path dataName u y t rate
end

%% AVERAGING (CONCATENATING)

% find dataObjs w/ matching expressions in titles and combine their data
inputExprs = ["A1","A2","E"];
titles = [dataObjs.title];

% initialize array of combined data objects
combinedObjs = [];

% for each input,
for idxInput = 1:length(inputExprs)
    % get objs containing relevant input
    mask = contains(titles,inputExprs(idxInput)); % idxs w/ this input
    relevantObjs = dataObjs(mask);

    % initialize new combined obj
    clear obj
    obj.u = [];
    obj.y = [];
    obj.t = [];

    % for each obj containing relevant input,
    for idxObj = 1:length(relevantObjs)
        relObj = relevantObjs(idxObj);

        % joining data
        obj.path(idxObj) = relObj.path;
        t_ = split(relObj.title,'_R'); % remove run # from title
        obj.title(idxObj) = t_(1);
        obj.u = [obj.u;relObj.u]; % concatenate
        obj.y = [obj.y;relObj.y]; % concatenate
        if(isempty(obj.t))
            obj.t = relObj.t;
        else
            obj.t = [obj.t;obj.t(end)+relObj.t]; % concatenate with addition
        end
        obj.rate(idxObj) = relObj.rate;
    end

    % add this combnied obj to the list
    combinedObjs = [combinedObjs,obj];
end

% replace previous array of dataObjs from files with new combined/concatenated ones
dataObjs = combinedObjs;
clear combinedObjs

%% PLOT TIME-SERIES
for idxObj = 1:length(dataObjs)
    t = dataObjs(idxObj).t;
    u = dataObjs(idxObj).u;
    y = dataObjs(idxObj).y;

    % initailize figure
    fig = figure;
    tl = tiledlayout(4,1);

    % labeling
    title(tl,dataObjs(idxObj).title(1),'Interpreter','none')
    xlabel(tl,'time (s)')
    
    % 4 inputs
    nexttile(1)
    plot(t,u(:,1:4))
    legend({'ail1 (deg)','ail2 (deg)','elev (deg)','gust (deg)'})
    title('Input')
    grid on
    set(gca,'XTickLabel',{})
    
    % 4 input derivatives
    nexttile(2)
    plot(t,u(:,5:8))
    legend({'ail1 (deg)','ail2 (deg)','elev (deg)','gust (deg)'})
    title('Input First Derivative')
    grid on
    set(gca,'XTickLabel',{})

    % 4 input double derivatives
    nexttile(3)
    plot(t,u(:,9:12))
    legend({'ail1 (deg)','ail2 (deg)','elev (deg)','gust (deg)'})
    title('Input Second Derivative')
    grid on
    set(gca,'XTickLabel',{})
    
    % outputs
    nexttile(4)
    plot(t,y)
    legend({'acc1 (g)','acc2 (g)','acc3 (g)','microstrain','pitch (deg)','pitch dot (deg/s)'})
    title('Output')
    grid on
end

%% COMPUTE FREQUENCY RESPONSE

% add path of czt_FRF() parent directory
addpath('../')

% czt_FRF() parameters
N = 10000;    % window size = 20 sec
w = [0.5,2]; % bandwidth
    % ^only go to 2 Hz because that's the highest frequency in the sweeps

for idxObj = 1:length(dataObjs)
    dataObj = dataObjs(idxObj);
    t = dataObj.t;
    for idxIn = 1:size(dataObj.u,2)
        u = dataObj.u(:,idxIn);
        for idxOut = 1:size(dataObj.y,2)
            y = dataObj.y(:,idxOut);
            % compute FRF components using czt_FRF()
            [~,dataObj.freq(:,idxIn,idxOut),...
                dataObj.coh(:,idxIn,idxOut),...
                Gxx_hat,Gyy_hat,Gxy_hat,Gyx_hat]...
                = czt_FRF(u,y,N,w,'off',dataObj.rate(1));

            % compute FRFs
            dataObj.H1_FRF(:,idxIn,idxOut) = Gxy_hat./Gxx_hat;
            dataObj.H2_FRF(:,idxIn,idxOut) = Gyy_hat./Gyx_hat;
        end
    end
    % "robust" Hv FRF
    dataObj.Hv_FRF = 0.5*(dataObj.H1_FRF+dataObj.H2_FRF);
    dataObjs_(idxObj) = dataObj;
end
dataObjs = dataObjs_;

%% SAVE FRF DATA
for idxObj = 1:length(dataObjs)
    dataObj = dataObjs(idxObj);
    savename = strcat('FRF_',strrep(dataObj.title(1),' ','_'));
    save(savename,'dataObj')
    disp(strcat('saved ',savename))
end

%% PLOT FRF DATA

% define I/O variable names
outputNames = ["acc1 (g)","acc2 (g)","acc3 (g)","microstrain","pitch (deg)","pitch dot (deg/s)"];
inputNames = ["ail1 (deg)","ail2 (deg)","elev (deg)","gust (deg)",...
    "ail1 dot (deg/s)","ail2 dot (deg/s)","elev dot (deg/s)","gust dot (deg/s)",...
    "ail1 ddot (deg/s2)","ail2 ddot (deg/s2)","elev ddot (deg/s2)","gust ddot (deg/s2)"];

% for each data object,
for idxObj = 1:length(dataObjs)
    dataObj = dataObjs(idxObj);
    
    % only use indices of nonzero control
    inputIdxs = find(logical(var(dataObj.u)));
    
    % define convenience variables
    nOutputs = size(dataObj.y,2);
    nInputs = size(dataObj.u,2);

    % initialize figure
    fig = figure;
    fig.Position = [0,0,1200,900];
    tl = tiledlayout(nOutputs,numel(inputIdxs));
    title(tl,dataObj.title(1),'FontSize',24,'Interpreter','none');
    xlabel(tl,'Frequency (Hz)','FontSize',18);
    ylabel(tl,'Magnitude (dB)','FontSize',18);
    
    % for each output,
    for idxOut = 1:nOutputs
        % for each input,
        for idxIn = inputIdxs
            freq = dataObj.freq(:,idxIn,idxOut);
            mag_H1 = abs(dataObj.H1_FRF(:,idxIn,idxOut));
            mag_H2 = abs(dataObj.H2_FRF(:,idxIn,idxOut));
            mag_Hv = abs(dataObj.Hv_FRF(:,idxIn,idxOut));

            % convert magnitude to dB
            mag_H1 = mag2db(mag_H1);
            mag_H2 = mag2db(mag_H2);
            mag_Hv = mag2db(mag_Hv);

            % initialize tile
            ax = nexttile();
            hold on
              % set(ax,'XScale','log')
            grid on
            
            % plot
            plot(freq,mag_Hv,'DisplayName','Hv FRF');
            % hold on
            plot(freq,mag_H1,'DisplayName','H1 FRF');
            plot(freq,mag_H2,'DisplayName','H2 FRF');

            % label I/O on top/left plots
            if(ax.Layout.Tile<=3)
                title(ax,inputNames(idxIn),'FontSize',12);
            end
            if(mod(ax.Layout.Tile,3)==1)
                ylabel(ax,outputNames(idxOut),'FontSize',12);
            end
        end
    end
end
legend()