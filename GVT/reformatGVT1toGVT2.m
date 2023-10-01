% Title: Re-format GVT1 Data to GVT2 Format
% Author: Anthony Su
% Date: 2023-07-27

close all
clear all
clc

% find GVT1 Data
dirFiles = dir("GVT1");
dirFileNames = strings(size(dirFiles));
for idxFile = 1:length(dirFiles)
    dirFileNames(idxFile) = dirFiles(idxFile).name;
end
GVT1FileNames = dirFileNames(contains(dirFileNames,'impact_') & ...
    contains(dirFileNames,'sens_') & ...
    contains(dirFileNames,'.mat'));
GVT1Paths = strings(size(GVT1FileNames));
for idxFile = 1:length(GVT1FileNames)
    GVT1Paths(idxFile) = string(['GVT1/',char(GVT1FileNames(idxFile))]);
end

% convert GVT1 Data
for idxPath = 1:length(GVT1Paths)
    % convert
    path = GVT1Paths(idxPath);
    newObj = GVT1toGVT2(path);

    % read I/O locations
    inputLoc = strings(0);
    outputLoc = strings(0);
    for idxSensor = 1:length(newObj.sensors)
        if(newObj.sensors(idxSensor).IsInput==true)
            inputLoc = [inputLoc,newObj.sensors(idxSensor).Location{:}];
        elseif(newObj.sensors(idxSensor).IsInput==false)
            outputLoc = [outputLoc,newObj.sensors(idxSensor).Location{:}];
        end
    end
    assert(length(inputLoc)<=1,'converted data has more than one input!')
    assert(length(inputLoc)>=1,'converted data has less than one input!')

    % new filename
    newName = ['TIME_impact_',char(inputLoc),'_sens_'];
    for idxOutput = 1:length(outputLoc)
        newName = [newName,char(outputLoc(idxOutput)),'-'];
    end
    newName(end) = []; % remove trailing hyphen
    newPath = ['GVT1_reformatted/',newName];

    % save
    rate = newObj.rate;
    sensors = newObj.sensors;
    H = newObj.H;
    save(newPath,"rate","sensors","H");
    disp(['saved ',newName])
end