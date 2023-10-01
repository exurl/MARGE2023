% Title: MARGE GVT
% Author: Anthony Su
% Date: 2023-05-20

close all
clear all

%% Initialize DAQ

% connect to DAQ
DAQ_avail = daqlist;
DAQ_obj = daq("ni");
DAQ_obj.Rate = 6400; % Hz
    % ^51200 Hz maximum

% import sensor specs and connect to sensors
sensorSpecs = readtable('sensor_data.xlsx');
for idx = 1:height(sensorSpecs)
    s.Name = sensorSpecs{idx,1};
    s.IsInput = sensorSpecs{idx,2};
    s.Unit = sensorSpecs{idx,3};
    s.Scaling= sensorSpecs{idx,4};
    s.Slot = sensorSpecs{idx,5};
    s.Channel = sensorSpecs{idx,6};
    if(s.IsInput)
        l_ = input(['where is the input "',s.Name{:},'" located?\n'],'s');
        s.Location = {l_};
    else
        l_ = input(['where is sensor "',s.Name{:},'" located?\n'],'s');
        s.Location = {l_};
    end
    sensors(idx) = s;
    ch(idx) = addinput(DAQ_obj,DAQ_avail{s.Slot,2},['ai',num2str(s.Channel)],'IEPE');
end

% check that there is only one input (hammer)
inputCt = 0;
for idx = 1:length(sensors)
    if(sensors(idx).IsInput)
        inputCt = inputCt+1;
        inputIdx = idx;
    end
end
if(inputCt==0)
    error('list of sensors does not contain input device!')
elseif(inputCt>1)
    error('list of sensors contains multiple input devices!')
end

%% Initialize Test

% prompt: how many tests to average
nTests = input('How many impacts in this dataset? ');

% prompt: how long after impact should record data for
fmin = input('What is the lowest frequency (Hz) behavior you are trying to capture? ');
assert(fmin>0,'frequency must be greater than zero!')
dtDesired = (1/fmin)*4;
disp(['Required signal length is ',num2str(dtDesired),' seconds after impact'])

% prompt: desired impact hammer force
forceDesired = input('What is the target impact hammer force (N)? ');

%% Perform Tests

% prompt: run name
% runName = input('Name of this run (file name): ','s');
runName = ['TIME_impact_',sensors(1).Location{:},'_sens_'];
for idx_sensor = 2:length(sensors)
    if(idx_sensor==length(sensors))
        runName = [runName,sensors(idx_sensor).Location{:}];
    else
        runName = [runName,sensors(idx_sensor).Location{:},'-'];
    end
end
[Y_,M_,D_] = ymd(datetime("today"));
runName = [runName,'_',num2str(Y_),'-',num2str(M_),'-',num2str(D_)];

% run tests
idxTests = 1;
while(idxTests<=nTests)
    % prompts: user that data collection begins
    input('Press ''Enter'' to begin data recording...','s');
    for ct = 3:-1:1
        fprintf([num2str(ct),'...\n'])
        pause(1)
    end
    beep
    disp('Recording data...')
    
    % record data
    start(DAQ_obj,'continuous')
    input('Press ''Enter'' to stop...\n','s')
    stop(DAQ_obj)

    % read data
    data = read(DAQ_obj,'all');

    % scale data
    for idx = 1:width(data)
        data{:,idx} = data{:,idx}.*1000*sensors(idx).Scaling;
    end

    % find impact
    for idx = 1:width(data)
        [pvs_,pis_] = findpeaks(data{:,idx});
        [pv_,pi_] = max(pvs_);
        peakVals(idx) = pv_;
        peakIdxs(idx) = pis_(pi_);
    end

    % check that there is dt seconds of data after impact
    dt = length(data{:,inputIdx}(peakIdxs(inputIdx):end))/DAQ_obj.Rate;
    if(dt<dtDesired)
        warning(['Did not take enough data! Took ',num2str(dt),' seconds when ',num2str(dtDesired),' seconds is required']);
        disp('Re-doing test...')
        continue
    end

    % check that the impact force is within limits
    impactForce = peakVals(inputIdx);
    if(abs(impactForce-forceDesired)>10)
        warning(['Impact force was ',num2str(impactForce),' N while desired force is ',num2str(forceDesired),' N'])
    end

    % plot impact
    f = figure();
    f.Position(1) = 0;
    f.Position(2) = 50;
    f.Position(3) = f.Position(3)*3;
    f.Position(4) = f.Position(4)*2;
    tl = tiledlayout(1,2);
    for idx = 1:length(sensors)
        labs_{idx} = [sensors(idx).Location{:},' (',sensors(idx).Unit{:},')'];
    end

    nexttile(1) % impact close-up
    stackedplot(data.Time,data{:,:},'DisplayLabels',labs_)
    xlabel('time')
    title(tl,runName,'Interpreter','none')
    xlim([data.Time(peakIdxs(inputIdx))-seconds(0.05),data.Time(peakIdxs(inputIdx))+seconds(0.15)])

    nexttile(2) % full data
    stackedplot(data.Time,data{:,:},'DisplayLabels',labs_)
    xlabel('time')

    % prompt: does diagnostic plot look indicate good data?
    isGood = input('Is this hit good? 1=YES 0=NO: ');
    if(isGood==1)
        H(idxTests).data = data;
        H(idxTests).peakIdxs = peakIdxs;
        idxTests = idxTests+1;
    end

    % close diagnostic plot
    try
        close(f)
    catch % if already closed by user
    end
end

% save data
doSave = input('Save data? 0=NO, 1=yes: ');
if(doSave~=0)
    filename = [runName,'.mat'];
    rate = DAQ_obj.Rate;
    save(filename,'H','sensors','rate')
    disp(['saved ',filename])
end