% Title: Anthony ASE LTI SS System Generation Driver
% Author: Anthony Su
% Date: 2023-08-07

clear all
% close all

%% INDEPENDENT VARIABLES/PARAMETERS

% NASTRAN data folder
nastranInputDir = 'ASEInputData/';

% filename to save state-space model as
filename = "ASE_SS.mat";

% NASTRAN structure data properties
NS = 15;
NC = 4;

% structure model parameters
ns = 15;
nc = 4; % [ail1, ail2, elev, vane], radians
    % ^NOTE: there is no way to disable controls

% decoupled structure damping ratios
zeta=zeros(ns,1);
zeta(2:9) = [0.028,0.042,0.030,0.112,0.031,0.018,0.022,0.032];
    % ^from GVT, see ./GVT/modalDataSummary.xlsx

% wild (but conservative) guess damping for modes of unknown damping
zeta(1) = 0.3;
zeta(10:ns) = 0.01;

% take only zeta corresponding to ns
zeta = zeta(1:ns);

% do we want to model gusts?
gusts = false;

% aero model properties
b = 0.2; % m

% aero model parameters
nLag = 4;
nLagG = 32;

% viscosity corrections
staticAeroCorrection = [1,1; 1,1];
flapCorrections = [1,1,1,1];
    % for [ail1, ail2, elev, vane]

% flight condition
rho = 1.293; % kg/m3
q = [60 100 164 207 281 343];

% physical constants
g = 9.80665; % m/s2

% accelerometer node IDs
accIds = [2223,2423,3008];

% strain gauge node ID
strainId = [2003];

% pitch sensor node ID
pitchId = [2003];

% ------------------------------------------------------------------------+ 
%% MODEL TUNING CHANGES (PRE-GENERATION)
% ------------------------------------------------------------------------+

% ==== MODEL SIZE ====
% gusts on?
gusts = false;
% slow input data
nastranInputDir = 'ASEInputData_slow/';
nLag = 0;
nLagG = 6;
% small model (2-DOF 0-LAG)
ns = 2;
nLag = 0;
zeta = zeta(1:ns);

% ==== NATURAL FREQUENCY REPLACEMENT ====
omegan = zeros(1,ns+nc); % Hz
% pitching frequency
% omegan(1) = 
% damping frequency
% omegan(2) = 

% ==== DAMPING RATIO SCALING ====
% pitching damping
zeta(1) = zeta(1)*100;
% bending damping
zeta(2) = zeta(2)*1;

% ==== AERODYNAMICS ====
% static aero corrections
staticAeroCorrection = ones(ns);
staticAeroCorrection(1,1) = 0.9; % CM_alpha
staticAeroCorrection(1,2) = 1; % CM_eta
staticAeroCorrection(2,1) = 0.5; % CL_alpha
staticAeroCorrection(2,2) = 1.5; % CL_eta
% control surface aero corrections
flapCorrections(1) = 0.6; % ail1
flapCorrections(2) = 0.6; % ail2
flapCorrections(3) = 0.6; % elev
flapCorrections(4) = 4; % vane

%% INTERMEDIATE VARIABLES

% stiffness matrix replacement
    % NOTE: omegan = sqrt(k/m) --> k = m*omegan^2, where m=1 in FEM
kNew = omegan.^2;

% flight conditions
u = sqrt(2*q/rho);
nSpeeds = length(u);

%% PLANT MATRICES GENERATION
plantObj = Anthony_ASE_SS_plant_generation(nastranInputDir,NS,NC,ns,nc,nLag,nLagG,gusts,u,rho,b,zeta,staticAeroCorrection,flapCorrections,kNew);
Ap = plantObj.A;
Bpc = plantObj.Bc;
if(gusts)
    BG = plantObj.BG;
else
    BG = [];
end

% speeds vector
u = plantObj.u;

% max valid frequency
omegaMax = plantObj.omegaMax;

%% OUTPUT MATRICES GENERATION
outputObj = Anthony_ASE_SS_output_generation(NS,NC,ns,nc,gusts,nSpeeds,Ap,Bpc,BG,accIds,strainId,pitchId);
Cp = outputObj.C;
Dpc = outputObj.Dc;
if(gusts)
    DG = outputObj.DG;
else
    DG = [];
end

%% RE-ORDER INPUT VECTOR VARIABLES

% rearrange order of inputs for [Bc], [Dc] from {u1 u2 ... u1Dot u2Dot
% ... u1DDot u2DDot ...} to instead become {u1 u1Dot u1DDot u2 ...}
TInputOrder = zeros(3*nc,3*nc);
for idxInput = 1:nc
    TInputOrder(idxInput,3*(idxInput-1)+1) = 1;
    TInputOrder(nc+idxInput,3*(idxInput-1)+2) = 1;
    TInputOrder(2*nc+idxInput,3*(idxInput-1)+3) = 1;
end
Bpc = pagemtimes(Bpc,TInputOrder);
Dpc = pagemtimes(Dpc,TInputOrder);

%% CONSTRUCT PLANT STATE-SPACE OBJECT
sysPlant = ss(Ap,Bpc,Cp,Dpc);

%% ACTUATOR MATRICES GENERATION AND IMPLEMENTATION

% compute actuator matrices
returnObj = Anthony_ASE_SS_actuator_generation();
Aa = returnObj.A;
Ba = returnObj.B;
Ca = returnObj.C;
Da = returnObj.D;
    % order: u = [u1 u1Dot u1DDot u2 ...]'

% construct actuator state-space object
sysAct = ss(Aa,Ba,Ca,Da);
sysAct(:,:,1:length(u)) = sysAct; % copy for each speed

% chain/append plant system to actuator system
sys = series(sysAct,sysPlant);
A = sys.A;
Bc = sys.B;
C = sys.C;
Dc = sys.D;

%% CHANGE INPUT UNITS/SIGNS TO MATCH EXPERIMENTAL SETUP/DATA

% input units: all must be converted from degrees to radians
Bc = Bc*deg2rad(1);
Dc = Dc*deg2rad(1);

% input signs: all must be flipped to change from experimental reference
% frame to modeling (NASTRAN) reference frame
Bc = -Bc;
Dc = -Dc;

%% CHANGE OUTPUT UNITS/SIGNS TO MATCH EXPERIMENTAL SETUP/DATA

% accelerometers must be converted from m/s^2 to gs
C(1:3,:,:) = C(1:3,:,:)/g;
Dc(1:3,:,:) = Dc(1:3,:,:)/g;

% strain gauge must be converted from strain to microstrain
C(4,:,:) = C(4,:,:)*1e6;
Dc(4,:,:) = Dc(4,:,:)*1e6;
% strain gauge signs must be flipped to change from modeling (NASTRAN) to experiment reference frame
C(4,:,:) = -C(4,:,:);
Dc(4,:,:) = -Dc(4,:,:);
% strain gauge calibration/correction to bring beam theory prediction in line with measured values
C(4,:,:) = 0.75*C(4,:,:);
Dc(4,:,:) = 0.75*Dc(4,:,:);

% pitch and pitch rate must be converted from radians to degrees
C(5:6,:,:) = C(5:6,:,:)*rad2deg(1);
Dc(5:6,:,:) = Dc(5:6,:,:)*rad2deg(1);
% pitch and pitch rate signs must be flipped to change from modeling (NASTRAN) to experiment reference frame
C(5:6,:,:) = -C(5:6,:,:);
Dc(5:6,:,:) = -Dc(5:6,:,:);

% ------------------------------------------------------------------------+
%% MODEL TUNING CHANGES (POST-GENERATION)
% ------------------------------------------------------------------------+

% John's CIFER correction
% A(ns+1,2,:) = 1300*A(ns+1,2,:); % effect: increase 1.45 hump in pitch out
% A(ns+1,2,:) = -5*A(ns+1,2,:); % effect: negligible

% flip gust input sign
Bc(:,4,:) = -Bc(:,4,:);
Dc(:,4,:) = -Dc(:,4,:);

% flip strain output sign
C(4,:,:) = -C(4,:,:);
Dc(4,:,:) = -Dc(4,:,:);

%% EXPORT MATRICES
if(gusts)
    save(filename,'A','Bc','BG','C','Dc','DG','omegaMax','u','NS','NC','ns','nc','nLag','nLagG','b','staticAeroCorrection','flapCorrections','rho','g','accIds','strainId','pitchId')
else
    save(filename,'A','Bc','C','Dc','omegaMax','u','NS','NC','ns','nc','nLag','nLagG','b','staticAeroCorrection','flapCorrections','rho','g','accIds','strainId','pitchId')
end
disp(['model generated and saved at ',char(filename)])

%% CALL RESPONSE ANALYSIS SCRIPTS

% compute FRFs
margeResponse

% plot FRFs against experiment
margeFreqExperiment