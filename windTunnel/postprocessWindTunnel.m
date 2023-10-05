% Title: Import Gust Vane Doublet Data
% Author: Anthony Su
% Date: 2023-08-30

% import, cleanup, and FFT of wind tunnel time-series data

% close all
clear all
clc

%% SELECT FILES

% define permutations of run names
speeds = ["060","100","164","207","281","343"];
    % trueSpeeds = [68,94,169,209,283,343]; % -John via Slack, 2023-09-06
inputs = ["A1S","A2S","ES","GD"];
runs = ["R1","R2","R3"];

baseDir = "./2023 summer/";
for idxSpeed = 1:length(speeds)
    speed = speeds(idxSpeed);

    % find all files for this speed
    fileObjs{idxSpeed} = dir(strcat(baseDir,'q',speed,'/*/output.mat'));

    % remove zeroing runs
    mask = false(size(fileObjs{idxSpeed}));
    for idxFile = 1:length(fileObjs{idxSpeed})
        fileObj = fileObjs{idxSpeed}(idxFile);
        if(contains(fileObj.folder,'zero'))
            mask(idxFile) = true;
        end
    end
    fileObjs{idxSpeed}(mask) = [];

    % get combined paths
    for idxFile = 1:length(fileObjs{idxSpeed})
        fileObj = fileObjs{idxSpeed}(idxFile);
        paths(idxSpeed,idxFile) = string(strcat(fileObj.folder,'\',fileObj.name));
    end
end

% reshape paths to be 1-dimensional
paths = reshape(paths,1,[]);

% remove "missing" elements, remnants from original array
mask = ismissing(paths);
paths(mask) = [];

% remove one-off runs
oneOffs = ["q060\A1D4_R1","q060\A1S5_R4","q100\ED3_R1"];
for idx1 = 1:length(oneOffs)
    oneOff = oneOffs(idx1);
    paths(contains(paths,oneOff)) = [];
end

% take only runs with A1S, A2S, ES, or GD (which all speeds have)
mask1 = contains(paths,'A1S');
mask2 = contains(paths,'A2S');
mask3 = contains(paths,'ES');
mask4 = contains(paths,'GD');
paths = paths(mask1 | mask2 | mask3 | mask4);

% sort paths into alphabetical order
paths = sort(paths);

% check files exist
if(any(~isfile(paths)))   
    error('file does not exist!')
end

% group paths
paths = reshape(paths,3,4,6);
    % paths(k,j,i) is speed i, input j, run k 

%% IMPORT, TRUCATE, REFORMAT
for idxSpeed = 1:length(speeds)
for idxInput = 1:length(inputs)
    % load data
    r1 = load(paths(1,idxInput,idxSpeed)).data;
    r2 = load(paths(2,idxInput,idxSpeed)).data;
    r3 = load(paths(3,idxInput,idxSpeed)).data;
    data = [r1,r2,r3];

    % data name
    pathParts = split(paths(1,idxInput,idxSpeed),["\","_"]);
    pathPart1 = pathParts(end-3); % e.g. q060
    pathPart2 = pathParts(end-2); % e.g. A1S5
    pathPart3 = pathParts(end-1); % e.g. R1

    dataName = pathPart2; % input, magnitude. e.g. A1S6
    
    temp = char(pathPart1);
    q = str2num(temp(2:4)); % speed
    
    temp = char(pathPart3);
    run = num2str(temp(2:end));

    for idxRun = 1:3
        %% IMPORT
        d = data(idxRun);
        
        % input data
        u_elev = get(d,'cmd').Values.elevator.Data;
        u_ail1 = get(d,'cmd').Values.aileron_1.Data;
        u_ail2 = get(d,'cmd').Values.aileron_2.Data;
        u_gust = get(d,'gust cmd').Values.Data;
        if(length(u_gust)==length(u_ail1)-1) % if gust data missing a point
            u_gust = [u_gust;0]; % assume the missing point is the last one
        end
        u = [u_ail1,u_ail2,u_elev,u_gust];

        % time vector
        t = get(d,'cmd').Values.elevator.Time;
        rate = (length(t)-1)/(t(end)-t(1));
        
        % output data 
        theta = get(d,'Plant:1').Values.theta.Data;
        temp = diff(theta)*rate;
        thetaDot = ([temp(1);temp]+[temp;temp(end)])/2; % numerical derivative
        acc1 = get(d,'Plant:1').Values.accel_1.Data;
        acc2 = get(d,'Plant:1').Values.accel_2.Data;
        acc3 = get(d,'Plant:1').Values.accel_3.Data;
        strain = get(d,'microstrain').Values.Data;
        y = [acc1,acc2,acc3,strain,theta,thetaDot];

        %% SHIFT AND TRUCATE
    
        % find time range of interest
        inputMask = find(sum(u~=0,2)); % indices with nonzero input
        idx1 = inputMask(1)-1;         % start 1 datapoint before input
        idx2 = inputMask(end)+5*rate;  % end 5 s after input
        if(idx2>length(y))
            idx2 = length(y);
        end
        if(idx1==0)
            idx1 = 1;
        end
        inputRange = idx1:idx2;

        % truncate data
        t = t(inputRange);
        u = u(inputRange,:);
        y = y(inputRange,:);
    
        % set new start time to zero
        t = t-t(1);
    
        % zero-mean everything
        % u = u-mean(u); % do not zero input ! -2023/10/02 John's advice
        y = y-mean(y);
            % NOTE: if data is truncated again in the future, it will have
            % to be centered again as well
        
        % store as object
        obj.path = paths(idxRun,idxInput,idxSpeed);
        obj.title = dataName;
        obj.fullTitle = string(['q',char(num2str(q)),'_',char(dataName)]);
        obj.q = q;
        obj.run = run;
        obj.u = u;
        obj.y = y;
        obj.t = t;
        obj.rate = rate;
        dataObjs(idxRun,idxInput,idxSpeed) = obj;
        clear u y t
    end
end
end

%% SAVE TIME-DOMAIN DATA
for idxSpeed = 1:length(speeds)
for idxInput = 1:length(inputs)
    data = dataObjs(:,idxInput,idxSpeed);

    % filter the accel data
    for idxRun = 1:3
        data(idxRun).y(:,1:3) = accelFilter(data(idxRun).y(:,1:3),rate);
    end

    % save time-domain data
    savename = strcat('./wtData/','TIME_',data.fullTitle);
    save(savename,'data')
    disp(['saved ',char(savename)])
end
end
clear data

%% CONCATENATE TIME-DOMAIN DATA
% ^combine data of 3 runs of each speed-input pair

for idxSpeed = 1:length(speeds)
for idxInput = 1:length(inputs)
    % initialize combined objec
    clear obj
    combinedObj.t = [];
    combinedObj.u = [];
    combinedObj.y = [];

    % join data
    for idxRun = 1:length(runs)
        % load single run
        obj = dataObjs(idxRun,idxInput,idxSpeed);

        % join single run to combined data
        combinedObj.u = [combinedObj.u;obj.u];
        combinedObj.y = [combinedObj.y;obj.y];
        if(isempty(combinedObj.t))
            combinedObj.t = obj.t;
        else
            combinedObj.t = [combinedObj.t;combinedObj.t(end)+obj.t]; % concatenate with addition
        end
        
    end

    % other obj properties
    combinedObj.q = obj.q;
    combinedObj.title = obj.title;
    combinedObj.fullTitle = obj.fullTitle;
    combinedObj.rate = obj.rate;

    % store
    combinedObjs(idxInput,idxSpeed) = combinedObj;
end
end

% this is the new dataObjs from now on
dataObjs = combinedObjs;
clear combinedObjs;

%% PLOT TIME-SERIES DATA
for idxObj = 1:numel(dataObjs)
    % plot
    plotTimeObj(dataObjs(idxObj))

    % save plot
    print(['TIME_',char(dataObjs(idxObj).fullTitle),'.png'],'-dpng','-r300')
end

%% COMPUTE FRFS

% add path of czt_FRF() parent directory
addpath('../')

% czt_FRF parameters
N = 2500; % window = 5 seconds (data rate is 500 Hz)
w = [0.4,2]; % bandwidth, Hz

% compute FRFs
for idxSpeed = 1:6
for idxInput = 1:4
% for idxSpeed = 1:length(speeds)
% for idxInput = 1:length(inputs)
    frfObjs(idxInput,idxSpeed) = computeObjFRF(dataObjs(idxInput,idxSpeed),N,w,rate);
end
end

% this is the new dataObjs from now on
dataObjs = frfObjs;
clear frfObjs

%% SAVE FREQUENCY-RESPONSE DATA
for idxSpeed = 1:length(speeds)
for idxInput = 1:length(inputs)
    data = dataObjs(idxInput,idxSpeed);

    % filter accel time-series data
    data.y(:,1:3) = accelFilter(data.y(:,1:3),rate);

    % save FRF+TIME objs
    savename = strcat('./wtData/','FRF_',data.fullTitle);
    save(savename,'data')
    disp(['saved ',char(savename)])
end
end
clear data

%% PLOT FRF DATA
for idxSpeed = 1:length(speeds)
    speed = idxSpeed;

    % plot it
    plotFrfObj(dataObjs(:,idxSpeed))

    % save plot
    print(['FRF_q',char(speeds(idxSpeed)),'.png'],'-dpng','-r300')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotTimeObj(dataObj)
% plots time-response dataObj

    t = dataObj.t;
    rate = (length(t)-1)/(t(end)-t(1)); 
    u = dataObj.u;
    y = dataObj.y;

    inputNames = {'ail1 (deg)','ail2 (deg)','elev (deg)','gust (deg)'};
    outputNames = {'acc1 (g)','acc2 (g)','acc3 (g)','hundreds microstrain','pitch (deg)'};
    % outputNames = {'acc1 (g)','acc2 (g)','acc3 (g)','hundreds microstrain','pitch (deg)','pitch dot (deg/s)'};
    
    % microstrain --> hundreds microstrain
    y(:,4) = y(:,4)*0.01;

    % filter accels
    for idxAcc = 1:3
        y(:,idxAcc) = accelFilter(y(:,idxAcc),rate);
    end

    % initailize figure
    fig = figure;
    fig.Position = [0,0,1500,750];
    tl = tiledlayout(2,1);

    % labeling
    title(tl,dataObj.fullTitle,'Interpreter','none')
    xlabel(tl,'time (s)')
    
    % inputs
    nexttile(1)
    plot(t,u(:,1:4))
    legend(inputNames)
    title('Input')
    grid on
    set(gca,'XTickLabel',{})
    
    % outputs
    nexttile(2)
    % plot(t,y)
    plot(t,y(:,1:5)) % ignore pitch rate
    legend(outputNames)
    title('Output')
    grid on
end

%%
function plotFrfObj(dataObjs)
% INPUTS:
    % dataObj : array of 4 dataObj, each corresponding to different inputs
    %           in order of ail1, ail2, elev, gust
    
    % define convenience variables
    nOutputs = size(dataObjs(1).y,2);
    nInputs = length(dataObjs);

    % initialize figure
    fig = figure;
    fig.Position = [0,0,1200,900];
    tl = tiledlayout(nOutputs,nInputs);
    
    % figure labeling
    title(tl,['q=',num2str(dataObjs(1).q)],'FontSize',24,'Interpreter','none');
    xlabel(tl,'Frequency (Hz)','FontSize',18);
    ylabel(tl,'Magnitude (dB)','FontSize',18);

    % for each output,
    for idxOut = 1:nOutputs
    % for each input,
    for idxIn = 1:nInputs
        dataObj = dataObjs(idxIn);

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
        
        inputNames = {'ail1 (deg)','ail2 (deg)','elev (deg)','gust (deg)'};
        outputNames = {'acc1 (g)','acc2 (g)','acc3 (g)','microstrain','pitch (deg)','pitch rate (deg/s)'};

        % plot
        plot(freq,mag_Hv,'DisplayName','Hv FRF');
        hold on
        plot(freq,mag_H1,'DisplayName','H1 FRF');
        plot(freq,mag_H2,'DisplayName','H2 FRF');

        % label I/O on top/left plots
        if(ax.Layout.Tile<=4) % top plots
            title(ax,dataObj.title,'FontSize',12);
        end
        if(mod(ax.Layout.Tile,4)==1) % left plots
            ylabel(ax,outputNames(idxOut),'FontSize',12);
        end
    end
    end

    % legend
    lg  = legend(); 
    lg.Layout.Tile = 'North'; % place to the right of the tiledlayout

    linkaxes(tl.Children,'x') % unify x-axes
end

%%
function frfObj = computeObjFRF(dataObj,N,w,rate)
% INPUTs
    % dataObj : time-domain dataObj
    % N       : window size
    % w       : bandwidth
% OUTPUT :
    % frfObj : superset of original obj, has new FRF properties

    frfObj = dataObj;
    for idxIn = 1:size(dataObj.u,2)
    
        for idxOut = 1:size(dataObj.y,2)

            % accelerometer issue management ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
            if(any(idxOut==[1,2,3]))
                u = dataObj.u(:,idxIn);
                y = dataObj.y(:,idxOut);

                % constrain accelerometer freq response bandwidth
                w_acc = [1,2]; % Hz
    
                % bandpass 0.6 to 20 Hz of time data
                % u = bandpass(u,[0.6,20],dataObj.rate); % commented out on 2023-10-02
                % y = bandpass(y,[0.6,20],dataObj.rate);

                % bandpass 0.8 to 25 Hz of time data
                u = accelFilter(u,rate);
                y = accelFilter(y,rate);

                % discard first 4 seconds of data
                % u = u(4*round(dataObj.rate):end); % NO!!!! DONT DO THIS!!
                % y = y(4*round(dataObj.rate):end);

                % re-center data
                u = u-mean(u);
                y = y-mean(y);

                % compute FRF components using czt_FRF()
                [~,frfObj.freq(:,idxIn,idxOut),...
                    frfObj.coh(:,idxIn,idxOut),...
                    Gxx_hat,Gyy_hat,Gxy_hat,Gyx_hat]...
                    = czt_FRF(u,y,N,w_acc,'off',dataObj.rate); % w/ w_acc
            else
                u = dataObj.u(:,idxIn);
                y = dataObj.y(:,idxOut);

                % compute FRF components using czt_FRF()
                [~,frfObj.freq(:,idxIn,idxOut),...
                    frfObj.coh(:,idxIn,idxOut),...
                    Gxx_hat,Gyy_hat,Gxy_hat,Gyx_hat]...
                    = czt_FRF(u,y,N,w,'off',dataObj.rate); % w/ w
            end

            % compute FRFs
            frfObj.H1_FRF(:,idxIn,idxOut) = Gxy_hat./Gxx_hat;
            frfObj.H2_FRF(:,idxIn,idxOut) = Gyy_hat./Gyx_hat;
        end
    end
    % "robust" Hv FRF
    frfObj.Hv_FRF = 0.5*(frfObj.H1_FRF+frfObj.H2_FRF);
end

%%
function xNew = accelFilter(x,rate)
% x : time-series signal
    w = [0.8,25];
    [b_,a_] = butter(3,w/rate,'bandpass'); % see John DM 2023-10-02
        % NOTE: B_ and A_ are coefficients for a discrete (z-domain) TF
    xNew = filter(b_,a_,x);
        % NOTE: filter takes coefficients for a discrete (z-domain) TF
end