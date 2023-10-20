% Title: Anthony ASE LTI SS Observer Generation
% Author: Anthony Su
% Date: 2023-08-07

function returnObj = Anthony_ASE_SS_output_generation(NS,NC,ns,nc,gusts,nSpeeds,A,Bc,BG,accIds,strainId,pitchId)
% INPUTS
    % NS       : # structural states in NASTRAN data
    % NC       : # rigid control states in NASTRAN data
    % ns       : # structural states in model output
    % nc       : # rigid control states in model output
    % gusts    : boolean, whether gusts are modeled
    % nSpeeds  : # flight conditions
    % A        : plant dynamics matrices
    % Bc       : control input dynamics matrices
    % BG       : gust input dynamics matrices
    % accIds   : IDs of NASTRAN nodes with z-accelerometers
    % strainId : ID of NASTRAN node with y-strain sensor
    % pitchId  : ID of NASTRAN node with y-rotation sensor
% OUTPUTS
    % C  : plant observation matrices
    % Dc : control input feedthrough matices
    % DG : gust input feedthrough matrices

%% INPUT PROCESSING
assert(numel(NS)==1)
assert(numel(NC)==1)
assert(numel(ns)==1 && ns<=NS)
assert(numel(nc)==1 && nc<=NC)
assert(numel(gusts)==1 && class(gusts)=="logical")
assert(numel(nSpeeds)==1)
assert(numel(strainId)==1)
assert(numel(pitchId)==1)

%% IMPORT DATA

% import NASTRAN mode shapes
PHI = load('ASEInputData\PHI_T.mat').PHI_T;

% decompose mode shapes into structural and control components
PHIs = PHI(:,1:ns);
PHIc = PHI(:,NS+1:NS+nc);

% import NASTRAN node locations
nodeIdsLocs = load('ASEInputData\GRID_ID_XYZ.mat').GRID_ID_XYZ;

% import NASTRAN modal loads
loadsPHI = load('ASEInputData\FM.mat').FM;

% take only loads corresponding to ns modes
loadsPHI = loadsPHI(:,1:ns);

%% ACCELEROMETER DEFINITION

% check accelerometer nodes are in list of nodes
assert(all(ismember(accIds,nodeIdsLocs(:,1))));

for idxAcc = 1:length(accIds)
    % get idx of accelerometer in node list
    idxNode = find(accIds(idxAcc)==nodeIdsLocs(:,1));

    % get modal z-displacement of accelerometers
    AccPHIz(idxAcc,:) = PHI((idxNode-1)*6+3,:); % z-disp is 3rd DOF of node
end

% structural component of z-displacement of accelerometers
AccPHIzs = AccPHIz(:,1:ns);

%% STRAIN GAUGE DEFINITION

% wing spar properties
sparThick = 0.0015875; % m
Izz = 2.541e-11; % m4
E = 6.89e10; % Pa

% get modal x-bending moment at root strain gauge
bendPHIx = loadsPHI(4,:); % x-bending moment is 4th load DOF, N-m
    % ^note: this is only for the 15 structural modes (does not have control modes)

% strain gauge calibration scaling
% strainMultiplier = 0.75; % original "calibration"
strainMultiplier = 1; % modified 2023-10-19

%% PITCH SENSOR DEFINITION

% get idx of pitch sensor in node list
idxNode = find(pitchId==nodeIdsLocs(:,1));

% get modal y-rotation of pitch sensor
pitchPHIy = PHI((idxNode-1)*6+5,:); % y-rotation is 5th DOF of node, rad

% structural component of y-rotation of pitch sensors
pitchPHIys = pitchPHIy(:,1:ns);

%% GENERATE {y} = [C]{x} + [Dc]{uc} + [DG]{uG}

% total number of states, outputs, inputs
nStates = size(A,1);
nInputs = size(Bc,2);
nOutputs = length(accIds)+length(strainId)+2*length(pitchId);

% [Tdisp] (Eq. 1.67)
Tdisp = zeros(ns,nStates);
Tdisp(1:ns,1:ns) = eye(ns);

% [Tvel] (Eq. 1.69)
Tvel = zeros(ns,nStates);
Tvel(1:ns,ns+1:2*ns) = eye(ns);

% initialize state-space matrices
C = zeros(nOutputs,nStates,nSpeeds);
Dc = zeros(nOutputs,nInputs,nSpeeds);
if(gusts)
    DG = zeros(nOutputs,3,nSpeeds);
end

for idxSpeed = 1:nSpeeds

    % accelerometer output (Eq. 1.71)
    CAcc = AccPHIzs*Tvel*A(:,:,idxSpeed); % output: m/s^2
    
    % accelerometer feedthrough
    DcAcc = AccPHIzs*Tvel*Bc(:,:,idxSpeed); % output: m/s^2
    if(gusts)
        DGAcc = AccPHIzs*Tvel*BG(:,:,idxSpeed); % output: m/s^2
    end

    % strain gauge output
    strainPHIy = -(bendPHIx*(-sparThick/2))/(Izz*E); % output: m/m
        % yStrain = -Mx*(z)/(I*E) where z = -t/2
    strainPHIy = strainMultiplier*strainPHIy;
    CStrain = strainPHIy*Tdisp;

    % strain gauge feedthrough
    DcStrain = zeros(1,nInputs);
    if(gusts)
        DGStrain = zeros(1,3);
    end

    % pitch sensor output
    CPitch = pitchPHIys*Tdisp; % output: rad
    CPitchDot = pitchPHIys*Tvel; % output: rad/s

    % pitch sensor feedthrough
    DcPitch = zeros(1,nInputs);
    DcPitchDot = zeros(1,nInputs);
    if(gusts)
        DGPitch = zeros(1,3);
        DGPitchDot = zeros(1,3);
    end

    % combine output matrices
    C(:,:,idxSpeed) = [CAcc;CStrain;CPitch;CPitchDot];

    % combine feedthrough matrices
    Dc(:,:,idxSpeed) = [DcAcc;DcStrain;DcPitch;DcPitchDot];
    if(gusts)
        DG(:,:,idxSpeed) = [DGAcc;DGStrain;DGPitch;DGPitchDot];
    end
end

returnObj.C = C;
returnObj.Dc = Dc;
if(gusts)
    returnObj.DG = DG;
end

end