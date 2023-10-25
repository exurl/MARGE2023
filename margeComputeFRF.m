% Title: Simulate MARGE Response
% Author: Anthony Su
% Date: 2023-08-01

% converted from margeResponse.m into a function on 2023-10-23

function dataObjs = margeComputeFRF(sys,q,omegaVec,saveName)
% INPUTS
    % sys      : array of state-space model objects at various flight conditions
    % q        : dynamic pressures (flight conditions)
    % omegaVec : frequency domain (Hz)
    % saveName : name of model; leave blank to not save FRFs
% OUTPUTS
    % dataObjs : array of time+freq data objects at various flight conditions

    % input validation
    if(~exist('saveName','var'))
        saveName = char([]);
    end
    
    % reference variables
    nSpeeds = length(q);
    nStates = size(sys.A,1);
    nInputs = size(sys.B,2);
    nOutputs = size(sys.C,1);
    
    % I/O names
    inputNames = ["aileron 1","aileron 2","elevator","gust vane"];
    outputNames = ["wing LE acc","wing TE acc","tail acc","root strain","pitch","pitch first derivative"];
    
    %% COMPUTE FREQUENCY RESPONSE
    
    % bode of all inputs
    for idxSpeed = 1:nSpeeds
        sys_ = sys(:,:,idxSpeed);
        [mag(:,:,:,idxSpeed),phase(:,:,:,idxSpeed),] = bode(sys_,omegaVec*2*pi);
            % indices of mag and phase are (idxOutput,idxInput,:,idxSpeed)
    end
    phase = deg2rad(phase);
    
    %% GENERATE RETURN OBJECT
    
    % convert mag+phase into complex FRF
    FRF = mag.*exp(phase*1j);
    
    % reformat matrices
    FRF = permute(FRF,[3,2,1,4]); % reorder indices to (:,idxOutput,idxInput,idxSpeed)
    
    % duplicate frequency vector for each I/O combo:
    freq = zeros(length(omegaVec),nInputs,nOutputs);
    for idxIn = 1:nInputs
        for idxOut = 1:nOutputs
            freq(:,idxIn,idxOut) = omegaVec;
        end
    end
    
    for idxSpeed = 1:nSpeeds
        FRF_ = FRF(:,:,:,idxSpeed);
    
        % create struct
        obj.path = [];
        obj.title = "State Space Model";
        obj.u = [];
        obj.y = [];
        obj.t = [];
        obj.rate = [];
        obj.freq = freq;
        obj.coh = [];
        obj.H1_FRF = FRF_;
        obj.H2_FRF = FRF_;
        obj.Hv_FRF = FRF_;
    
        % store in array
        dataObjs(idxSpeed) = obj;
    end
    
    % save array of FRF data structs along with corresponding speeds
    if(~isempty(saveName))
        fileName = ['FRF_',char(saveName),'.mat'];
        save(fileName,'dataObjs')
        disp(['saved ',fileName]);
    end

end