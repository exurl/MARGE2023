% Title: Plot State-Space and Experimental FRFs Together
% Author: Anthony Su
% Date: 2023-08-17

% converted from margeFreqExperiment.m into a function on 10-24-2023

function [magError,phaseError] = margeCompareFRF(dataObjs)
% INPUTS
    % dataObjs : array of time+freq data objects at various flight conditions
% OUTPUTS
    % residual : a measure of error between state-space model and experiment
    % error    : error of each omega/y/u/q combo

    nSpeeds = 6;
    nInputs = 4;
    nOutputs = 6;
    nPoints = length(dataObjs(1).freq(:,1,1));

    %% IMPORT EXPERIMENTAL DATA
    % do this only once when called repeatedly

    persistent expObjs badData
    if(isempty(expObjs))
        % q60
        expObjs = load('windTunnel\wtData\FRF_q60_A1S5.mat').data;
        expObjs(2,1) = load('windTunnel\wtData\FRF_q60_A2S5.mat').data;
        expObjs(3,1) = load('windTunnel\wtData\FRF_q60_ES2.mat').data;
        expObjs(4,1) = load('windTunnel\wtData\FRF_q60_GD4.mat').data;
        
        % q100
        expObjs(1,2) = load('windTunnel\wtData\FRF_q100_A1S5.mat').data;
        expObjs(2,2) = load('windTunnel\wtData\FRF_q100_A2S5.mat').data;
        expObjs(3,2) = load('windTunnel\wtData\FRF_q100_ES2.mat').data;
        expObjs(4,2) = load('windTunnel\wtData\FRF_q100_GD4.mat').data;
        
        % q164
        expObjs(1,3) = load('windTunnel\wtData\FRF_q164_A1S5.mat').data;
        expObjs(2,3) = load('windTunnel\wtData\FRF_q164_A2S5.mat').data;
        expObjs(3,3) = load('windTunnel\wtData\FRF_q164_ES2.mat').data;
        expObjs(4,3) = load('windTunnel\wtData\FRF_q164_GD4.mat').data;
        
        % q207
        expObjs(1,4) = load('windTunnel\wtData\FRF_q207_A1S5.mat').data;
        expObjs(2,4) = load('windTunnel\wtData\FRF_q207_A2S5.mat').data;
        expObjs(3,4) = load('windTunnel\wtData\FRF_q207_ES2.mat').data;
        expObjs(4,4) = load('windTunnel\wtData\FRF_q207_GD4.mat').data;
        
        % q281
        expObjs(1,5) = load('windTunnel\wtData\FRF_q281_A1S5.mat').data;
        expObjs(2,5) = load('windTunnel\wtData\FRF_q281_A2S5.mat').data;
        expObjs(3,5) = load('windTunnel\wtData\FRF_q281_ES2.mat').data;
        expObjs(4,5) = load('windTunnel\wtData\FRF_q281_GD4.mat').data;
        
        % q343
        expObjs(1,6) = load('windTunnel\wtData\FRF_q343_A1S3p5.mat').data;
        expObjs(2,6) = load('windTunnel\wtData\FRF_q343_A2S3p5.mat').data;
        expObjs(3,6) = load('windTunnel\wtData\FRF_q343_ES1.mat').data;
        expObjs(4,6) = load('windTunnel\wtData\FRF_q343_GD4.mat').data;
    
        % take only relevant FRFs from each experiment
        for idxSpeed = 1:nSpeeds
            for idxInput = 1:nInputs
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

        % combine all objects
        for idxSpeed = 1:nSpeeds
        for idxInput = 1:nInputs
        for idxOutput = 1:nOutputs
            newExpObjs(idxSpeed).freq(:,idxInput,idxOutput) = expObjs(idxInput,idxSpeed).freq(:,idxOutput);
            newExpObjs(idxSpeed).H1_FRF(:,idxInput,idxOutput) = expObjs(idxInput,idxSpeed).H1_FRF(:,idxOutput);
            newExpObjs(idxSpeed).H2_FRF(:,idxInput,idxOutput) = expObjs(idxInput,idxSpeed).H2_FRF(:,idxOutput);
            newExpObjs(idxSpeed).Hv_FRF(:,idxInput,idxOutput) = expObjs(idxInput,idxSpeed).Hv_FRF(:,idxOutput);
            newExpObjs(idxSpeed).title(idxInput,idxOutput) = expObjs(idxInput,idxSpeed).title;
        end
        end
        end
        expObjs = newExpObjs;
        clear newExpObjs;

        %% MASK OF BAD DATA

        % initialize
        badData = false(6,4,6); % idxs are (output,input,speed)

        % ignore all accelerometers
        badData(1:3,:,:) = true;
            % ^NOTE: in the future when you remove this, make sure to add
            % specific bad accelerometer FRFs below

        % ignore pitchDot (since pitch is already accounted for)
        badData(6,:,:) = true;

        % q60 gust-->pitch
        badData(5:6,4,1) = true;
        
    end

    %% COMPUTE RESIDUAL

    % initialize error matrix
    magError = zeros(nPoints,nOutputs,nInputs,nSpeeds);
    phaseError = zeros(nPoints,nOutputs,nInputs,nSpeeds);

    % compute error
    for idxSpeed = 1:nSpeeds
    for idxInput = 1:nInputs
    for idxOutput = 1:nOutputs
    if(badData(idxOutput,idxInput,idxSpeed)==false)
        % extract values
        freqExp = expObjs(idxSpeed).freq(:,idxInput,idxOutput);
        frfExp = expObjs(idxSpeed).Hv_FRF(:,idxInput,idxOutput);
        freq = dataObjs(idxSpeed).freq(:,idxInput,idxOutput);
        frfSS = dataObjs(idxSpeed).Hv_FRF(:,idxInput,idxOutput);

        % interpolate experiment to get corresponding value
        frfExp = interp1(freqExp,frfExp,freq,'makima','extrap'); % if runtime is long try switching to linear interpolation

        % error
        magError(:,idxOutput,idxInput,idxSpeed) = ((abs(frfSS)-abs(frfExp))./abs(frfExp)).^2; % normalized magnitude error squared
        phaseError(:,idxOutput,idxInput,idxSpeed) = (phase(frfSS)-phase(frfExp)).^2; % phase error squared
    end
    end
    end
    end

    % phase error modulo 2pi
    phaseError = mod(phaseError,2*pi);
    phaseError(phaseError>pi) = 2*pi-phaseError(phaseError>pi);

end