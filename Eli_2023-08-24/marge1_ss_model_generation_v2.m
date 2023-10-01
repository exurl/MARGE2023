%
% State Space Model - LTI - of an aeroelastic plant 
% subject to control surface inputs (by irreversible controls) 
% and gust inputs
%

clear all
close all
clc
%%
% Inputs:
%
% NS - number of input structural modes
NS=15;%15
ns = 2;%number of output strutural modes
%
% NC - number of input rigid body control surface modes
NC=4;
nc = 4;%number of output structural modes
%
% nLag - number of lag terms in the Roger approximation 
% of the structural/control modes
nLag=0;%4

%
% nLagG - number of lag terms in the Roger approximation of the gust vector
%
nLagG=32;
nLagG = 0;
%
% b - refereNCe semichord
b=0.4/2;
b = 0.08318;%meters reference semichord
%
% Reduced FrequeNCies: kbar
load ASEInputData/k_bar.mat; kbar=k_bar;
%
% kMax  - the highest reduced frequeNCy used for Roger RFA fitting
kMax=max(kbar);
%
% NK - the size of reduced frequeNCies
NK=length(kbar);

% Pbar0, Pbar1, Pbar2, Pbar3, ... stored in a single array Pbar
% of dimeNSioNS (NS+NC) x (NS+NC) x (2+nLag)
% That is, Pbar ( (NS+NC)x(NS+NC), 1) is Pbar0,
%          Pbar ( (NS+NC)x(NS+NC), 2) is Pbar1, etc.
% For the case of irreversible controls we use only the (NS x (NS+NC) part
% of the matrices
load ASEInputData/A0.mat; load ASEInputData/A1.mat; load ASEInputData/A2.mat;
load ASEInputData/D.mat; load ASEInputData/E.mat; %load ASEInputData/R.mat

Pbar(:,:,1)=A0;
Pbar(:,:,2)=A1;
Pbar(:,:,3)=A2;
for i=1:nLag
	Pbar(:,:,i+3)=D(:,(NS+NC)*(i-1)+1:(NS+NC)*i)*E((NS+NC)*(i-1)+1:(NS+NC)*i,:);
end
%
% PbarG0, PbarG1, PbarG2, PbarG3, ... stored in a single array PbarG
% of dimeNSioNS (NS+NC) x (1) x (2+nLagG)
% That is, PbarG ( (NS+NC)x 1, 1) is PbarG0,
%          PbarG ( (NS+NC)x 1, 2) is PbarG1, etc.
% For the case of irreversible controls we use only the (NS x 1) part
% of the matrices
load ASEInputData/AG0.mat; load ASEInputData/AG1.mat; load ASEInputData/AG2.mat;
load ASEInputData/DG.mat; load ASEInputData/EG.mat; %load ASEInputData/RG.mat

PbarG(:,:,1)=AG0;
PbarG(:,:,2)=AG1;
PbarG(:,:,3)=AG2;
for i=1:nLagG
	PbarG(:,:,i+3)=DG(:,i)*EG(i,:);
end
%
% betabar - vector of lag terms used in the Roger RFA fit for structural
%           and control modes. DimeNSion nLag x 1
for i=1:nLag
    betabar(i)=1.7*kMax*i^2/(nLag+1)^2;
end
%
% betabarG - vector of lag terms used in the Roger RFA fit for the gust
%            vector. DimeNSion nLagG x 1
for i=1:nLagG
    betabarG(i)=1.7*kMax*i^2/(nLagG+1)^2;
end
%
% M - generalized mass matrix (NS+NC) x (NS+NC)
load ASEInputData/MHH_T.mat; M=MHH_T;
%
% K - generalized stiffness matrix (NS+NC) x (NS+NC)
load ASEInputData/KHH_T.mat; K=KHH_T;
%
% Eigenvactors: PHI_T N x (NS+NC) (each PHI_i is traNSformed into a vector {x1 y1 z1 r1 p1 q1 x2 y2 z2 r2 p2 q2 ...}')
load ASEInputData/PHI_T.mat; PHI=PHI_T;
%
% zeta - vector of structural damping ratios for the structural modes. 
%        DimeNSion : NS x 1
zeta = 0.01*ones(NS,1);
zeta(1) = 0.01;
zeta(2) = 0.0282;%measured 1st wing bending damping
%

% NrespoNSes - number of respoNSes to be studied
%
% for each respoNSe - a vector ModalContributioNS iNCluding the modal 
%                     contributioNS to the respoNSe.
%                     For example, if the respoNSe is a vertical
%                     displacement at a structural node, then the vector
%                     ModalContributioNS would contain the displacements of
%                     the point in each of the mode shapes

%% Additional Matrices

% ss, cc, and sc matrices
M_ss=M(1:ns,1:ns);
M_cc=M(NS+1:end,NS+1:end);
M_sc=M(1:ns,NS+1:end);

K_ss=K(1:ns,1:ns);
K_cc=K(NS+1:end,NS+1:end);
K_sc=K(1:ns,NS+1:end);
PHI = PHI(:,[1:ns (NS+1):end]);
for i=1:nLag+3
	Pbar_ss(:,:,i)=Pbar(1:ns,1:ns,i);
	Pbar_cc(:,:,i)=Pbar(NS+1:end,NS+1:end,i);
	Pbar_sc(:,:,i)=Pbar(1:ns,NS+1:end,i);
end

for i=1:nLagG+3
	PbarG_s(:,:,i)=PbarG(1:NS,1,i);
	PbarG_c(:,:,i)=PbarG(NS+1:end,1,i);
end

% Create a structural damping matrix (approximated) based on the diagonal
% terms of [M] and [K]
% Note: when the natural modes of the structure are used K_ss and M_ss will
% be diagonal. But the formulation here allows for the usage of mode shapes
% that are not natural, and in that case M_ss and K_ss would not be diagonal.
%
% Note for a SDOF oscillator, the equation is writtn:
%
% m xdotdot(t) + c xdot(t) + k x(t) = F(t)
%
% xdotdot(t) + 2.*zeta*wn*xdot(t)+ wn^2*x(t) = F(t)/m
% where wn is the natural frequency.
% With this notation: c = m*2.*zeta*wn
%
% Find the uncoupled natural frequencies:
%
for i=1:ns
  OmegaUncoupled(i)=sqrt(K_ss(i,i)/M_ss(i,i));
end
%
% print uncoupled natural frequencies, if desired
%
' Uncoupled Structural Frequencies (Hz) ';
OmegaUncoupledHz=OmegaUncoupled/2./pi;
OmegaUncoupledHz
% 
% Create the structural C matrix as a diagonal matrix made of
% 2 * zeta(i) * OmegaUncoupled(i) * M_ss(i,i) terms 
% 
%
% Note: The formulation at this stage is limited to a diagonal damping
% matrix even when M_ss and K_ss are not diagonal.
%
Css=zeros(ns,ns);
%
for i=1:ns
  Css(i,i)=2.*zeta(i,1)*OmegaUncoupled(i)*M_ss(i,i);
end
%% Extract aerodynamic and gust matrices if needed

% Aerodynamic Matrices: Q
load ASEInputData/QHH_T.mat; Q=QHH_T;

% Gust Vectors: QG
load ASEInputData/QHG_T.mat; QG=QHG_T;

for i=1:length(kbar)
	Q_ss(:,:,i)=Q(1:ns,1:ns,i);
	Q_cc(:,:,i)=Q(NS+1:end,1:NS+1:end,i);
	Q_sc(:,:,i)=Q(1:ns,1:NS+1:end,i);

	QG_s(:,:,i)=QG(1:ns,1,i);
	QG_c(:,:,i)=QG(NS+1:end,1,i);
end

%% Accelerometers (ACC)
% Loading Grid ID and their locatioNS
load ASEInputData/GRID_ID_XYZ.mat;

% ACC Grid ID, the order is forward accel, aft accel, tail accel
GRID_ID_ACC=[2223,2623,3008];

% Finding locatioNS ACC grid points
N_ACC=length(GRID_ID_ACC);
for i=1:N_ACC
    ind_ACC(i)=find(GRID_ID_ACC(i)==GRID_ID_XYZ);
end
%ind_ACC is the node number. Given the node number, we need to find the
%corresponding row in the modal matrix. We want the T3 (z-displacement) row
%and we want the T5 (pitch rotation) row.

moment_arm = [.08 -.08 .67];%rotation moment arm for acc1, acc2, acc3

% Finding Z-displacement for each mode
for i=1:N_ACC
	
    T1_row_index = (ind_ACC(i)-1)*6 + 1;%each node has 6 DOF (x,y,z,roll,pitch,yaw)
    T3_row_index = T1_row_index + 2;%1 + 2 = 3
    T5_row_index = T1_row_index + 4;%1 + 4 = 5
    T3_row = PHI(T3_row_index,1:ns);
    T5_row = PHI(T5_row_index,1:ns);
    PHI_Z_ACC(i,1:ns) = T3_row + moment_arm(i)*T5_row;%PHI_Z_ACC(i,1:ns+NC) = T3_row + moment_arm(i)*T5_row;

end
%% PITCH/ANGLE-OF-ATTACH DOF
% Loading Grid ID and their locatioNS
load ASEInputData/GRID_ID_XYZ.mat;

% ACC Grid ID, the order is forward accel, after accel, tail accel
GRID_ID_PITCH=[2003];%root location

% Finding locatioNS pitch grid points
N_PITCH=length(GRID_ID_PITCH);
for i=1:N_PITCH
    ind_PITCH(i)=find(GRID_ID_PITCH(i)==GRID_ID_XYZ);
end
%ind_PITCH is the node number. Given the node number, we need to find the
%corresponding row in the modal matrix. We want the T5 (pitch rotation) row.


% Finding pitch rotation for the node at the root
for i=1:N_PITCH
	
    T1_row_index = (ind_PITCH(i)-1)*6 + 1;%each node has 6 DOF (x,y,z,roll,pitch,yaw)
    
    T5_row_index = T1_row_index + 4;%1 + 4 = 5
    
    T5_row = PHI(T5_row_index,1:ns);%T5_row = PHI(T5_row_index,1:ns + NC);
    PHI_PITCH(i,1:ns) = T5_row;%PHI_PITCH(i,1:NS+NC) = T5_row;

end











%
% ---------------------
% - Flight ConditioNS -
% ---------------------
%
% Sea Level air deNSity 
%

%Air density used for a 90F day
AirDeNSity=1.15 ;% kg/m^3 %     < -------------------------------AirDeNSity
rho=AirDeNSity;
%
% State space models can be generated at a single speed or multiple speeds.
% For flutter analysis we need a fine enough array of speeds.
%
% For flutter analysis (by another code) NSpeeds speeds are used 
% to span a range of speeds from SpeedLow to SpeedHigh. 
% At each speed a state space model is
% generated, to be used later for flutter analysis by finding the
% eigenvalues of the [Ap] matrix as fuNCtioNS of speed.
%                         <--------------------One Speed -------NSpeeds=
%
%  SpeedUnitConversionFactor - used to convert speeds range input, given in
%  user selected units, to the units coNSistent with the units of the
%  structural dynamic and aerodynamic inputs.
%
% For speed in m/s:
%
SpeedUnitConversionFactor=1. ; % < -------------------------------------------UunitConversion
%
% Find speeds from dynamic pressures in N/m^2
%

%vector of dynamic pressures
q_vector = [68 94 169 209 283 343];
NSpeeds = length(q_vector);

DynPressLow = min(q_vector); % < --------------------------------------------------------------qDLow

DynPressHigh = max(q_vector); % < -------------------------------------------------------------qDHigh


% speedlow = 10;%m/s
% 
% 
% speedHigh = 24;%m/s




% Speed step size (for the case of multiple speeds)
% Here: Dynamic Pressure step size
if(NSpeeds > 1)
%dspeed=(speedHigh-speedLow)/(NSpeeds-1);
dDynPressure=(DynPressHigh-DynPressLow)/(NSpeeds-1);
else
%dspeed=0;
dDynPressure=0;
end
%
% ---------------------------------------------------------------
%
% Prepare state space model matrices (eq. 1.38): [Ap], [Bc], {BG}
% at either a single flight speed or a set of speeds.
%
% Allocate storage space for arrays
%
% Total number of states
%
NStatesTotal=2*ns+ns*nLag+nLagG;
NStates=NStatesTotal;
%
% an array for storing the speeds used
%
Speeds=zeros(NSpeeds,1);
DynamicPressures=zeros(NSpeeds,1);
%
% Eqs. (1.38)
%
% an array storing the [Ap] matrices for the different speeds
ApMatrices=zeros(NStates,NStates,NSpeeds);
% an array storing the [Bc] matrices for the different speeds
BcMatrices=[];
if(NC > 0);
BcMatrices=zeros(NStates,3*NC,NSpeeds);
end
% an array storing the [BG] matrices for the different speeds
BGMatrices=[];
igust = 'no';%no gusts
if(igust == "YesGusts");
BGMatrices=zeros(NStates,3,NSpeeds);
end
%
% Array for storing values of highest frequenCies (rad/sec) for which the
% math models are valid
%
OmegaMax=zeros(NSpeeds,1);
%
% initialize Uspeed before beginning the loop over speeds
%
%Uspeed=speedLow-dspeed;
DynPress=DynPressLow-dDynPressure;
%
for ispeed=1:NSpeeds                                %     <------ Loop on Dyn Pressures
    %
    % Uspeed is in units coNSistent with the units used for analysis.
    % SpeedUnitConversionFactor convert speeds to units desired for
    % presentation of results
    %
    %Uspeed=Uspeed+dspeed;
    %*******************************
    %DynPress=DynPress+dDynPressure;

    DynPress = q_vector(ispeed);

    DynamicPressures(ispeed)=DynPress;
    %
    Uspeed=sqrt(DynPress*2./rho);
    %
    % Speed values are stored in their coNSistent units
    Speeds(ispeed)=Uspeed;
    % Uflight=Uspeed used alternatively. To be "cleaned" later for more
    % compact programming.
    Uflight=Uspeed;
    %
    % Dynamic Pressure
    %
    %qD(ispeed)=0.5*rho*Uspeed^2;
    qD(ispeed)=DynPress;
    %
%
% Find the highest frequNCy (rad / sec) for which the model is valid:
%
OmegaMax(ispeed)=kMax*Uspeed/b;
%
%
% create Mbarbarss,Cbarbarss, and Kbarbarss (equatioNS 1.17, 1.18)
% only the NS x NS part is needed
%
Mbarbarss=zeros(ns,ns);
Cbarbarss=zeros(ns,ns);
Kbarbarss=zeros(ns,ns);
%note that:
% Pbar_ss(:,:,1) = P0 or A0 stiffness
% Pbar_ss(:,:,2) = P1 or A1 damping
% Pbar_ss(:,:,3) = P2 or A2 mass
%
Mbarbarss(:,:)=M_ss(:,:)-0.5*rho*b*b*Pbar_ss(:,:,3);%good
Kbarbarss(:,:)=K_ss(:,:)-0.5*rho*Uflight*Uflight*Pbar_ss(:,:,1);%good
Cbarbarss(:,:)=Css(:,:)-0.5*rho*Uflight*b*Pbar_ss(:,:,2);%good
%
Mbarbarsc=[];
Cbarbarsc=[];
Kbarbarsc=[];
%
if (NC > 0)
%
% create Mbarbarsc,Cbarbarsc, and Kbarbarsc (equatioNS 1.17, 1.18)
% only the NS x NC part is needed
%
%
Mbarbarsc=zeros(ns,nc);
Cbarbarsc=zeros(ns,nc);
Kbarbarsc=zeros(ns,nc);
%
Mbarbarsc(:,:)=M_sc(:,:)-0.5*rho*b*b*Pbar_sc(:,:,3);%good
Kbarbarsc(:,:)=-0.5*rho*Uflight*Uflight*Pbar_sc(:,:,1);%good
Cbarbarsc(:,:)=-0.5*rho*Uflight*b*Pbar_sc(:,:,2);%good
% Note: no [Ksc] and [Csc] here.
%
%
% end if loop on whether there are control surfaces (if NC ne 0)
else
end
%
% Find the inverse of Mbarbarss (required to convert eq. 1.28 to standard
% state-space form)
%
InvMbarbarss=zeros(ns,ns);
%
InvMbarbarss=inv(Mbarbarss);
%
%
% Convert the Lag terms from their Roger fit values 
% to their Laplace domain values (eq. 1.12)
%
if(nLag >0)
for i=1:nLag;
  beta(i)=betabar(i)*Uspeed/b;
end
%
 ' Aero Lag Roots - beta '
for ii=1:nLag
    LagRoots(ii,1)=beta(ii);
end
LagRoots
else
end
%
%
betaG=[];
if(igust == "YesGusts")
for i=1:nLagG;
  betaG(i)=betabarG(i)*Uspeed/b;
end
%
 ' Aero Gust Lag Roots - betaG '
for ii=1:nLagG
    GustLagRoot(ii,1)=betaG(ii);
end
GustLagRoot
% 
else
end
%
% Create the aerodynamic lag equatioNS for the structural and control
% motioNS (eq. 1.23, 1.24)
%
Ar=[];
Brs=[];
Brc=[];
if(nLag > 0)
    %
Ar=zeros(ns*nLag,ns*nLag);
Brs=zeros(ns*nLag,ns);
if(NC > 0);
Brc=zeros(ns*nLag,nc);
else
end
%
for iLags=1:nLag;
  istart=(iLags-1)*ns+1;
  iend=iLags*ns;
  %
  for i=istart:iend
  Ar(i,i)=-beta(iLags);
  end
%
  for i=1:ns
  for j=1:ns
Brs(istart-1+i,j)=Pbar_ss(i,j,3+iLags);
  end    
  end
%
if (nc > 0)
  for i=1:ns
  for j=1:nc
  Brc(istart-1+i,j)=Pbar_sc(i,j,3+iLags);
  end    
  end
else
end
%
% end loop on iLags
end
%
else
end
%
%
ArG=[];
BrG=[];
%
if (igust == "YesGusts")
%
% Create the aerodynamic lag equatioNS for the gust column (eqs. 1.25, 1.26)
%
if(nLagG > 0)
ArG=zeros(nLagG,nLagG);
BrG=zeros(nLagG,1);
%
  for i=1,nLagG;
  ArG(i,i)=-betaG(i);
  BrG(i,1)=1.;  
  end
%
else
end
else
end
%
% Prepare the matrix [Ir} (eqs. 1.28 and 1.29)
%
Ir=[];
if( nLag > 0)
%
Ir=zeros(ns,ns*nLag);
for iLags=1:nLag;
jstart=(iLags-1)*ns;
for i=1:ns
Ir(i,jstart+i)=1.;
end  
end
else
end
%
%
PGr=[];
%
if (igust == "YesGusts")
% Prepare the matrix PGr (eqs. 1.28,1.30)
%
if(nLagG >0 )
PGr=zeros(ns,nLagG);
  for i=1:ns
  for j=1:nLagG
    PGr(i,j)=PimbarGs(i,1,j+3);
  end
  end
else
end
else
end
%
% Eqs. 1.32
%
T21=-InvMbarbarss*Kbarbarss;
T22=-InvMbarbarss*Cbarbarss;
%
%**********************************************************************
T2r=[];
if(nLag >0 )

T2r=0.5*rho*Uspeed*Uspeed*InvMbarbarss*Ir; % < ------------------------------ 
else
end
% *********************************************************************
%
%***********************************************************************
T2Gr=[];
if(igust == "YesGusts")
    if(nLag >0 )
T2Gr=0.5*rho*Uspeed*InvMbarbarss*PGr;
    else
    end
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
T2c=[];
%
if (NC > 0)
%
% T2c
%
Temp=zeros(ns,3*nc);
%
  for i=1:ns
  for j=1:nc
  Temp(i,j)=Kbarbarsc(i,j);    
  Temp(i,j+nc)=Cbarbarsc(i,j);
  Temp(i,j+2*nc)=Mbarbarsc(i,j);
  end
  end
%
T2c=zeros(ns,3*nc);
T2c=-InvMbarbarss*Temp;
%
clear Temp;
%
else
end
%
%
T2G=[];
%
if (igust == "YesGusts")
    %
% Prepare T2G
%
Temp=zeros(ns,3);
%
Temp(:,1)=PbarGs(:,1,1); % PbarGs0
Temp(:,2)=PbarGs(:,1,2); % PbarGs1
Temp(:,3)=PbarGs(:,1,3); % PbarGs2
%
T2G=zeros(ns,3);
%
T2G=0.5*rho*Uspeed*InvMbarbarss*Temp;  %%%%%%% Multiply by b? ????????????????
%
clear Temp;
%
else
end
%
% < ---------------------------------------------------------------single speed model
%
% Prepare the plant state space matrices at a single speed
% to be stored in the arrays of plant matrices at all speeds
% This code is derived from a code capable of generating state space
% models for multiple velocities
% 
% Eqs. 1.34-1.38
%
% Total number of states = NStates
%
Ap=zeros(NStates,NStates);
Bc=[];
if (nc >0)
Bc=zeros(NStates,3*nc);
else
end
BG=[];
if(igust == "YesGusts")
BG=zeros(NStates,3);
else
end
%
% The Ap matrix
% -------------
%
for i=1:ns
Ap(i,ns+i)=1.;
end
%
for i=1:ns
  for j=1:ns
Ap(i+ns,j)=T21(i,j);
Ap(i+ns,j+ns)=T22(i,j);    
  end
end
%
if( nLag > 0)
for i=1:ns
  for j=1:ns*nLag
Ap(i+ns,2*ns+j)=T2r(i,j);    
  end
end
else
end
%
%
if (igust == "YesGusts") 
for i=1:ns
  for j=1:nLagG
Ap(i+ns,2*ns+nLag*ns+j)=T2Gr(i,j);    
  end
end
else
end
%
if(nLag > 0) 
for i=1:nLag*ns %
  for j=1:ns
Ap(2*ns+i,ns+j)=Brs(i,j);    
  end
end
%
for i=1:nLag*ns
  for j=1:nLag*ns
Ap(2*ns+i,2*ns+j)=Ar(i,j);    
  end
end
else
end
%
if (igust == "YesGusts")
    if(nLagG > 0)
for i=1:nLagG;
  for j=1:nLagG;
 Ap(2*ns+nLag*ns+i,2*ns+nLag*ns+j)=ArG(i,j);   
  end
end
    else
    end
else
end
%
% End of Ap matrix generation
%
%
if (nc > 0)
%    
% Generate [Bc]
% -------------
%
Bc=zeros(NStates,3*nc);
%
for i=1:ns
  for j=1:3*nc
    Bc(ns+i,j)=T2c(i,j);
  end
end
%
if(nLag >0)
for i=1:nLag*ns
  for j=1:ns
    Bc(2*ns+i,nc+j)=Brc(i,j);
  end
end
else
end
%
% End of [Bc] generation
%
else
end
%
%
if(igust == "YesGusts")
%    
% Generate [BG]
% -------------
%
BG=zeros(NStates,3);
%
for i=1:ns
  for j=1:3
BG(NS+i,j)=T2G(i,j);    
  end
end
%
if(nLagG > 0)
for i=1:nLagG
BG(2*ns+nLag*ns+i,2)=BrG(i,1);  
end
else
end
%
% End of [BG] generation
%
else
end
%
% The ASE model - the state equatioNS:
%
% {xdot]=[Ap]{x}+[Bc]{uc}+[BG]{uG}
%
% {uc} =   {qc}
%          s{qc}
%         s^2{qc}
%
% {uG} = wG
%       s wG
%      s^2 wG
%
% Store in the arrays that contain system matrices for all speeds
% The Ap,Bc,BG matrices for the current speed ispeed are stored in the
% arrays ApMatrices, BcMatrices, BGMatrices, where the third index identify
% the speed case.
%
ApMatrices(:,:,ispeed)=Ap;
%
if(nc > 0)
BcMatrices(:,:,ispeed)=Bc(:,:);    
else
end
%
if(igust == "YesGusts")
BGMatrices(:,:,ispeed)=BG(:,:);
else
end
%
% end loop over speeds
%
end
%
clear Ap Bc BG ;
%
% End of generation of Ap, Bc, and BG matrices.

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%OUTPUT EQUATIONS SHOULD GO BELOW THIS LINE
%
% Convert acceleration outputs to g's
%
% gravitational acceleration in m/s^2
gravitation = 9.81; % < ---------------------------------------gravity
%                                                                          Note Units
PHIsACC(1:N_ACC,1:ns)=PHI_Z_ACC(:,1:ns);
% Provide accelerations in g's
PHIsACC=PHIsACC/gravitation;

Tdisp=zeros(ns,NStates);
Tvel=zeros(ns,NStates);
%
for ii=1:ns
Tdisp(ii,ii)=1.;
Tvel(ii,ii+ns)=1.;
end

% Loop over speeds again
%
for ispeed=1:NSpeeds
%
% {Vector of accelerations}=[Cp]{xp}+[Dcp]{uc}+DGp]{ug}
%
Ap(:,:)=ApMatrices(:,:,ispeed);
%
CpACC=PHIsACC*Tvel*Ap;
%
DcpACC=[];
Bc=[];
if(nc > 0)
Bc(:,:)=BcMatrices(:,:,ispeed);    
DcpACC=PHIsACC*Tvel*Bc;
end
%
DGpACC=[];
BG=[];
if(igust == "YesGusts")
BG(:,:)=BGMatrices(:,:,ispeed);
DGpACC=PHIsACC*Tvel*BG;
end
% Vector of displacements
CpDisplacements = PHI_PITCH(:,1:ns)*Tdisp;


%======================== STRAIN OUTPUT EQUATION ========================
% Vector of root loads
%
% Loads at the root FM(1:6,i) all 6 components for a mode i. 
% Z force for all modes is FM(3,:) 
load C:\Users\John\Documents\Code\Eli_code\Eli_input_data\FM.mat; % <-------------------------------------------------FM
%                                                                    ---------Root
% {RootLoads}=[FM}{qs}=[FM][Tdisp]{xp}
%
FMs=FM(:,1:ns);
CpRootLoads=FMs*Tdisp;
%FM gives x,y,z forces and x,y,z moments
%we want Mx, the moment about the x-axis. (which runs along the fuselage).
%A positive bending moment should bend the wing CW when viewed from the
%front using this convention. This may or may not be the correct
%convenction for marge.
CpRootLoads = CpRootLoads(4,:);

%stress = My/I and stress = E*strain, thus strain = stress/E = My/(EI). For
%marge, y = thickness/2, E = 6.78e10, I = 0.2499 in^4, thickness = t = 1/16
%inch
E_strain = 6.78e10;%Pa
I_strain = 0.2499*(0.0254)^4;%m^4
t_strain = 0.0015875;%meter, thickness
moment_to_strain_conversion = (t_strain/2)/(E_strain*I_strain);
CpRootLoads = CpRootLoads*moment_to_strain_conversion;
nLoads=0;
%============================END STRAIN OUTPUT EQUATION ================

Cp = [CpDisplacements;CpACC;CpRootLoads];


if(NC >0)
Dp = [zeros(size(CpDisplacements,1),3*nc);DcpACC;zeros(size(CpRootLoads,1),3*nc)];
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%???
%
BpcMatrices(:,:,ispeed)=Bc(:,:);
CpMatrices(:,:,ispeed)=Cp(:,:);
sizeC=size(CpMatrices);
DpcMatrices(:,:,ispeed)=Dp(:,:);
DpGMatrices(:,:,ispeed)=DGpACC(:,:);
%
% end loop over speeds
%
end

%
%
% Ref area = 12 ft^2, span=85 iNChes, mac = ???
%
Sref=12.*144.; % < -----------------------------------------------------------Sref
span=85.; % < ----------------------------------------------------------------span
%span;
mac=1.; % < ------------------------------------------------------------------mac
%

% nLOADSoutputs=nLoads;
nLOADSoutputs = 999;%fake number to make the code run ++++++++++++++++++++++++++++++++++++++++++++++++++
%
NConditioNS=NSpeeds;
%
save('marge1_initial_models.mat', ...
     'NS', 'NC', 'NStates','nLOADSoutputs','NConditioNS', ...
     'Speeds', 'DynamicPressures','OmegaMax', ...
     'ApMatrices', 'BpcMatrices', 'CpMatrices','DpcMatrices',...
     'AirDeNSity','Sref','span','b');
                         
%% plotting
r_start = 0;
g_start = 255;
b = 0;
size_start = 50;%10
size_end = 10;
for jj = 1:NSpeeds
    r = r_start + jj*8;
    g = g_start - jj*8;
    b = 0;
    color = [r g b]/255;
    subplot(1,2,1)
    scatter(real(eig(ApMatrices(:,:,jj))),imag(eig(ApMatrices(:,:,jj))),NSpeeds - (jj-1),color,'filled')
    hold on
    grid on
    xlabel('Imaginary Part')
    ylabel('Real Part')
    ylim([0 80])
    xlim([-10 10])

    subplot(1,2,2)
    scatter(real(eig(ApMatrices(:,:,jj))),imag(eig(ApMatrices(:,:,jj))),NSpeeds - (jj-1),color,'filled')
    hold on
    grid on
    xlabel('Imaginary Part')
    ylabel('Real Part')
    ylim([2 10])
    xlim([-1.2 0.4])
end
% xlabel('Imaginary Part')
% ylabel('Real Part')
% ylim([0 80])
% xlim([-10 10])

