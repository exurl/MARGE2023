%
% State Space Model - LTI - of an aeroelastic plant
% subject to control surface inputs (by irreversible controls)
% and gust inputs
%
%
%
% No Lags
%
% originally written by Eli Livne, modified by Anthony Su
%
% April 13 2023 - Model for Aeroelastic Response to Gust Excitation (MARGE)
%
% <TODO>: double check corrections
%
% 3 accelerometers, 1 root load, 1 hall effect pitch sensor
%
% Generate a state space model at a number of conditions, to be used later
% for control law synthesis
%
% First the aeroelastic state space model (no actuators and sensors)
%
% The ASE LTI State Space Model model:
% ------------------------------------
%
% {xdot]=[Ap]{x}+[Bc]{uc}+[BG]{uG}
%
% Control inputs in the form of forced control surface rotations and their
% first and second time derivatives (Eq. 1.27,1.34):
%
% {uc} =   {qc}
%          s{qc}
%         s^2{qc}
%
% Gust inputs in the form of the amplitude of an upward-directed vertical
% gust velocity as a function of time and its first and second time
% derivatives (Eq. 1.27,1.34).
%
% {uG} = wG
%       s wG
%      s^2 wG
%
%
% Notes on numerical inputs:
%
%
% The code is written for length units in m, force units in N and time in s
%
%
% Changing to other units requires changes in internal unit conversion factors
% that would take user-defined inputs and converrt them to consistemt
% units for simulation.
%
% **************
% * Inputs:    *
% **************
%
% Numbers of structural modes and rigid-body control modes
% ---------------------------------------------------------
%
% NS - number of input structural modes (from NASTRAN. The structural modes
% are for the airframe with the control surfaces locked, connected to the
% airframe, via spring that represent the stiffnesses of locked actuators.
%
% NC - number of input rigid body control surface modes (from NASTRAN). The
% "rigid" control rotation modes describe rotations of the control surfaces
% one by one on their hinges as rigid bodies, free to rotate. These motions
% will be driven by the control system and actuators in the closed loop
% case.
%
% In the MARGE & gust generation case,
% one control mode is the rotation of the gust vanes, which is treated as
% just another rotation (to be commanded) of a control surface.
%
% Unsteady Aerodynamics
% ----------------------
%
% nLag - number of lag terms in the Roger approximation
%        of the aerodynamic matrix terms associated with
%        structural/control modes
%
% nLagG - number of lag terms in the Roger approximation of the gust
%         force vector
%
% Geometry for unsteady aerodynamic scaling:
%
% b - reference semichord
%
% kMax  - The highest reduced frequency used for
%         Roger (Rational Function Approximation) RFA fitting (eq. 1.2)
%
% Roger Approximation Pbar and PbarG matrices:
%
% Pbar0, Pbar1, Pbar2, Pbar3, ... (eq. 1.2) stored in a single array Pbar
% of dimensions (NS+NC) x (NS+NC) x (3+nLag)
% That is, Pbar ( (NS+NC)x(NS+NC), 1) is Pbar0,
%          Pbar ( (NS+NC)x(NS+NC), 2) is Pbar1, etc.
% For the case of irreversible controls we use only the (NS x (NS+NC) part
% of the matrices
%
% The code is written currently for one gust input. The aerodynamic gust
% matrices are, therefore, column vectors.
%
% PbarG0, PbarG1, PbarG2, PbarG3, ...(eq. 1.13)stored in a single array PbarG
% of dimensions (NS+NC) x (1) x (3+nLagG)
% That is, Pbar ( (NS+NC)x 1, 1) is Pbar0,
%          Pbar ( (NS+NC)x 1, 2) is Pbar1, etc.
% For the case of irreversible controls we use only the (NS x 1) part
% of the matrices
%
% betabar - vector of lag terms used in the Roger RFA fit for structural
%           and control modes. Dimension nLag x 1 (eq. 1.2)
%
% betabarG - vector of lag terms used in the Roger RFA fit for the gust
%            vector. Dimension nLagG x 1 (eq. 1.13). Note the shift from
%            un-barred and barred elements of the Roger Gust Force
%            approximation is the same as the transformation used for the
%            ss (structure-structure) aerodynamic force terms.

%
% Flight conditions:
% ------------------
%
% rho - Atmospheric density
% Uspeed - speed at which the state space model is created
%           This can be a single speed (if a single state space model is
%           desired at some given speed) or a set of speeds (for flutter
%           analysis).
%
% Ulow, Uhigh: range of speeds for flutter root locus analysis, if
% required.
% If a state space model at a single speed is desired, then enter
% Ulow=Uhigh=The Desired Speed.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
clear all
close all
%
format long g;
%
%   Inputs
% ---------------------------
%
% Geometry (inch) - the reference semi-chord used for the reduced frequency
%
b=0.4/2; %     < ------------------------------------------------------------b
%
% NS - number of input structural modes (from NASTRAN)
%
NS=15;%        < -------------------------------NASTRAN DATA ------------------NS
%
% NC - number of input rigid body control surface modes (from NASTRAN)
%
NC=4;%        < --------------------------------NASTRAN DATA ------------------NC
%
% In the case here there are 4 rigid body modes: 2 rigid rotations of the
% wing flaps on their hinges, one rigid rotation of the elevator
% and one rigid rotation of the gust vanes
% (that move together).
%
% The number of modes actually used (a subset of all available modes):
%
% Less than NS structural modes can be used.
% The user may experiment with models based
% on increasing numbers of modes used.
%
%
ns=NS; %------------Used to generate the LTI model-----------------------------ns
%
%
% The number of control surfaces actually used out of the NC available
% nc
%
nc=NC; % < ---------Used to generate the LTI model-------------------------------nc
%
%
% The order of appearance of control modes in the control modes set is
% defined later. There are 2 TE flap modes, One elevator mode and one
% gust vanes modes.
%
%
%
%
%
%
% ---------------------------
% - Structural Dynamics     -
% ---------------------------
%
% M - Generalized mass matrix (NS+NC) x (NS+NC) from NASTRAN
% K - Generalized stiffness matrix (NS+NC) x (NS+NC) from NASTRAN
%     The partitions of the K matrix associated with rigid control modes
%     should be zero
%
%
%
load ../ASEInputData/MHH_T.mat; %     < ------------------------------------MHH
M=MHH_T;
%
% K - generalized stiffness matrix (NS+NC) x (NS+NC)
%
load ../ASEInputData/KHH_T.mat; %     < ------------------------------------KHH
K=KHH_T;
%
% zeta - vector of structural damping ratios for the structural modes.
%        Dimension : NS x 1
% zeta is a viscous damping coefficient. The equivalent structural
% hysteretic damping coefficient would be g = 2 * zeta
%
zeta=zeros(NS,1);
zeta(2:9) = [0.028,0.042,0.030,0.112,0.031,0.018,0.022,0.032]; % <---------- zeta
    % ^from GVT, see ./GVT/modalDataSummary.xlsx

%
% The ss and sc mass and stiffness matrices (eqs. 1.14 and 1.16)
% Those are extracted from the NASTRAN generated matrices.
%
Mss=M(1:ns,1:ns);
%
Mcc=M(NS+1:NS+nc,NS+1:NS+nc);
Msc=M(1:ns,NS+1:NS+nc);
Mcs=Msc';

Kss=K(1:ns,1:ns);
% We apply a correction (actually, replace the NASTRAN generate value)
% of Kss(1,1) to give us the <TODO> Hz frequency for the 1st bending mode
% This is with a generalized mass matrix that is the unit matrix.
% This needs to change when many cases are run, to work with the actual
% NASTRAN based M and K matrices.
%
%
%
% Because rigid-body rotations of control surfaces for control encounter no
% structural stiffness (the control surfaces are free to rotate
% structurally, and it would be the actuators that keep them in place) the
% cs and cc parts of the stiffness matrix are not used.
% They should, actually, be zero.
%
% Find the coupled natural frequencies with the control surfaces locked in
% place on their actuator springs (The "structural" modes)
% That is, the static stiffness of actuators is represented by rotation
% springs.
%
% [Kss]{eigvector}=lambda[Mss]{eigenvector}, with lambda=omega^2
%
lambda=eig(Kss,Mss);
lambda=sqrt(lambda);
% convert frequencies from rad/sec to Hz
lambda=lambda/2./pi; % NOTE TODO BUG: this should be divided by 4*pi^2
% sort and order frequencies from low to high (in Hz)
frequencies=sort(lambda);
% print natural frequencies
" Natural Frequencies (Hz)"
frequencies
%
% Convert highest natural frequency back to (rad/sec)
% It will be used to find the highest reduced frequency
%
wHighest=max(lambda)*2*pi;
%
clear lambda
%
% Create a structural damping matrix (approximated) based on the diagonal
% terms of [M] and [K]
% Note: when the natural modes of the structure are used Kss and Mss will
% be diagonal. But the formulation here allows for the usage of mode shapes
% that are not natural, and in that case Mss and Kss would not be diagonal.
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
for i=1:ns;
    OmegaUncoupled(i)=sqrt(Kss(i,i)/Mss(i,i));
end
%
% print uncoupled natural frequencies, if desired
%
' Uncoupled Structural Frequencies (Hz) ';
OmegaUncoupledHz=OmegaUncoupled/2./pi;
OmegaUncoupledHz
%
% Create the structural C matrix as a diagonal matrix made of
% 2 * zeta(i) * OmegaUncoupled(i) * Mss(i,i) terms
%
%
% Note: The formulation at this stage is limited to a diagonal damping
% matrix even when Mss and Kss are not diagonal.
%
Css=zeros(ns,ns);
%
for i=1:ns;
    Css(i,i)=2.*zeta(i,1)*OmegaUncoupled(i)*Mss(i,i);
end
%

% ------------------
% Aerodynamics     -
% ------------------
%
% nLag - number of lag terms in the Roger approximation
% of the structural/control modes (Eqs. 1.2 and 1.11)
%
nLag=4; %     < --------------------------------------------------------------nLag
%
% Are gust inputs included? Yes: igust='YesGusts', No: igust='NoGusts'
%
% nLagG is the number of aerodynamic lag states for the unsteady gust
% forces.
% we start with nLagG=0 in case the case does not include gust inputs
% This will be overridden in cases where gust inputs are included.
%
nLagG=0; %
%
igust="YesGusts"; %     < -----------------------------------------------------igust
%
%
if(igust == "YesGusts")
    %
    % nLagG - number of lag terms in the Roger approximation of the gust vector
    % (used only when gust vector Roger matrices are provided
    % and if igust = 'YesGusts' (Eq. 1.13)
    %
    % Note: often a large number of lag terms is required to cature the
    % frequency dependency of the gust force vector accurately by a Roger RFA.
    %
    nLagG=32; %     < ------------------------------------------------------------nLagG
else
end
%
%
% Tabulated Reduced Frequencies: kbar (reduced frequencies for which
% unsteady aero matrices were available from NASTRAN, to be fit by Roger
% RFA.
% The Roger RFA fitting is done before this code is used and its matrices
% are provided as inputs in a directory: ASEInputData
%
% The tabulated reduced frequencies are not used here, except for finding
% kMax - the highest reduced frequency used for the Roger fit.
%
% kbar and k_bar store the tabulated reduced frequencies.
%
load ../ASEInputData/k_bar.mat;  %     < ---------------------------------k_bar
kbar=k_bar;
%
% kMax  - the highest reduced frequency used for Roger RFA fitting
%
%< ---------------------------------------------------------------------------kMax
kMax=max(kbar)
%
% The highest reduced frequency used for the Roger RFA generation tells us
% what the highest actual frequency is (omega) above which the state space
% model is not valid.
%
% < ------------------------ a note about the potential addition (later) of a low-pass filter
% Since k = omega * b /U, then omega = k * U /b (rad/sec)
% and the max valid omega: omegaMax = kMax * U / b
% The state space model should not be used above this frequency at the
% speed U for which it is generated. A low-pass filter may be added to the
% system to filter out behavior above omegaMax.
%
% NK - the number of tabulated reduced frequencies (generated by NASTRAN).
%
NK=length(kbar);
%
% Pbar0, Pbar1, Pbar2, Pbar3,(eq. 1.2)... (eqs. 1.2, 1.12)
% stored in a single array Pbar
% of dimensions (NS+NC) x (NS+NC) x (3+nLag)
% That is, Pbar ( (NS+NC)x(NS+NC), 1) is Pbar0,
%          Pbar ( (NS+NC)x(NS+NC), 2) is Pbar1, etc.
% For the case of ireversible controls we use only the (NS x (NS+NC) part
% of the matrices
load ../ASEInputData/A0.mat; %     < --------------------------------------------A0
load ../ASEInputData/A1.mat; %     < --------------------------------------------A1
load ../ASEInputData/A2.mat; %     < --------------------------------------------A2
D=[];
E=[];
if(nLag >0)
    load ../ASEInputData/D.mat;  %     < --------------------------------------------D
    load ../ASEInputData/E.mat;  %     < --------------------------------------------E
else
end
% load ../ASEInputData/R.mat  % not used.
% [R] is used when using Minimum State RFA (Karpel's Method).
%
% Pbar is of full dimension (NS+NC) x (NS+NC) x (3+nLag)
% by "full dimension" we mean the data provided by NASTRAN
%
Pbar(:,:,1)=A0(:,:); % Pbar0
Pbar(:,:,2)=A1(:,:); % Pbar1
Pbar(:,:,3)=A2(:,:); % Pbar2
%
% The matrices used for the analysis here are sub-matrices of the full
% matrices (for just the number of modes and controls we are interested
% in).
%
% The ss (structure-structure) aero Roger matrices (Eq. 1.16)
%
Pbarss(:,:,1)=A0(1:ns,1:ns);
Pbarss(:,:,2)=A1(1:ns,1:ns);
Pbarss(:,:,3)=A2(1:ns,1:ns);
%
% Pbar(3) to Pbar(3+nLag): The lag matrices
% from the Nastran model by Marat Mor (only if there are lag terms)
if (nLag > 0)
    for i=1:nLag
        Pbar(:,:,i+3)=...
            D(:,(NS+NC)*(i-1)+1:(NS+NC)*i)*E((NS+NC)*(i-1)+1:(NS+NC)*i,:);
    end
    %
    % The ss structure-structure lag matrices
    % These are (ns x ns x (3+nLag)
    %
    for i=1:nLag
        % ss terms of the lag matrices
        Pbarss(:,:,i+3)=Pbar(1:ns,1:ns,i+3);
    end
    %
else
end
%
% Corection factor for all aerodynamic boxes <
% < --------------------------------------------------------------------------AeroCorrection
%
% Thickness Effect on cLalpha
AeroCorrection=1.00; % < ------------------------------------------------------0.92 <TODO>: adjust these
%
Pbarss=Pbarss*AeroCorrection ;
%
% Check if there are control surfaces. If so extract the sc aero matrices
% eqs. 1.14, 1.16
% They are of dimension ns x NC x (3+nLag)
% Remember that NSused=ns
%
Pbarsc=[];
FlapCorrections=[];
%
if(nc > 0)
    %
    Pbarsc(:,:,1)=A0(1:ns,NS+1:NS+nc);
    Pbarsc(:,:,2)=A1(1:ns,NS+1:NS+nc);
    Pbarsc(:,:,3)=A2(1:ns,NS+1:NS+nc);
    %
    if(nLag > 0) ;
        for i=1:nLag
            % sc terms of the lag matrices
            Pbarsc(:,:,i+3)=Pbar(1:ns,NS+1:NS+nc,i+3);
        end
    else
    end
    %
    % Correction factors for control surface aerodynamics:
    %
    % Correct {qc} <
    % ---------------------------------------------------------------------------{qc}
    % ---------------------------------------------------------------------------correction
    FlapCorrections=zeros(nc,nc); % < -------------------------------------------enter
    %
    %
    %
    %
    % there used to be LARGE's control surfaces here
    %
    %
    %
    FlapCorrections(1,1) = 1.0; % onboard flap 0.7  <TODO>: correct
    FlapCorrections(2,2) = 1.0; % outboard flap 0.7 <TODO>: correct
    FlapCorrections(3,3) = 1.0; % root rigid body rotation
    FlapCorrections(4,4) = 1.0; % gust vanes as control surfaces
    %
    %
    for ii=1:3+nLag
        %
        % The following can be done with a single Matlab command instead of the
        % double loop
        %    Temp(:,:)=Pbarsc(:,:,ii);
        for i=1:ns
            for j=1:nc
                Temp(i,j)=Pbarsc(i,j,ii);
            end
        end
        Temp=Temp*FlapCorrections;
        Temp=Temp*AeroCorrection;
        % AeroCorrection needed for Qsc on top of control corrections? < ----------
        %
        % Enter the corrected aerodynamic matrix control columns into Pbarsc
        %    Pbarsc(:,:,ii)=Temp;
        for i=1:ns
            for j=1:nc
                Pbarsc(i,j,ii)=Temp(i,j);
            end
        end
    end
    %
else
end
%
clear Temp;
%
% Gust input
%
PbarG=[];
if (igust == "YesGusts")
    %
    % PbarG0, PbarG1, PbarG2, PbarG3, ... stored in a single array PbarG
    % of dimensions (NS+NC) x (1) x (3+nLagG)
    % That is, PbarG ( (NS+NC)x 1, 1) is PbarG0,
    %          PbarG ( (NS+NC)x 1, 2) is PbarG1, etc.
    % For the case of irreversible controls we use only the (NS x 1) part
    % of the matrices. That is: a single gust velocity input wG.
    load ../ASEInputData/AG0.mat; %     < -------------------------------------------AG0
    load ../ASEInputData/AG1.mat; %     < -------------------------------------------AG1
    load ../ASEInputData/AG2.mat; %     < -------------------------------------------AG2
    load ../ASEInputData/DG.mat;  %     < -------------------------------------------DG
    load ../ASEInputData/EG.mat;  %     < -------------------------------------------EG
    %load ../ASEInputData_NoLags/RG.mat
    %
    % PbarG: NS x 1 x (3+nLagG)
    %
    PbarG(:,1,1)=AG0(1:(NS+NC),1);
    PbarG(:,1,2)=AG1(1:(NS+NC),1);
    PbarG(:,1,3)=AG2(1:(NS+NC),1);
    %
    if (nLagG > 0)
        %
        for i=1:nLagG
            PbarG(:,:,i+3)=DG(:,i)*EG(i,:);
        end
    else
    end
    %
    % The structure component of PbarG
    %
    % PbarGs: ns x 1 x (3+nLagG)
    %
    PbarGs(:,1,1)=AG0(1:ns,1);
    PbarGs(:,1,2)=AG1(1:ns,1);
    PbarGs(:,1,3)=AG2(1:ns,1);
    %
    if (nLagG > 0)
        %
        for i=1:nLagG
            % The lag vectors of PbarG (eq. 1.13)
            PbarGs(:,1,i+3)=PbarG(1:ns,1,i+3);
        end
    else
    end
    %
    % Correct PbarGs
    %
    % < --------------------------------------------------------------------------Aero Correction
    % < --------------------------------------------------------------------------Gust columns
    PbarGs=PbarGs*AeroCorrection ;
    %
    % end if loop that checks whether gust loads are included
else
end
%
% betabar - vector of lag roots used in the Roger RFA fit for structural
%           and control modes. Dimension nLag x 1
% Eqs. 1.2, 1.5
%
if (nLag > 0)
    %
    for i=1:nLag
        betabar(i)=1.7*kMax*i^2/(nLag+1)^2; %     < ------------------------------betabar
    end
    %
    ' betabar '
    betabar
else
end
%
betabarG=[];
if (igust == 'YesGusts')
    %
    % betabarG - vector of lag roots used in the Roger RFA fit for the gust
    %            vector. Dimension nLagG x 1
    % Eq. 1.13
    %
    if (nLagG > 0)
        for i=1:nLagG
            betabarG(i)=1.7*kMax*i^2/(nLagG+1)^2; %     < ----------------------------betabarG
        end
        %
        ' betabarG '
        betabarG
    else
    end
    %
else
end
%
% ---------------------
% - Flight Conditions -
% ---------------------
%
% Sea Level air density in lbf, inch, sec units (snail/inch^3)
%
AirDensity=1.293 ; % in kg/m^3 %     < -------------------------------AirDensity
rho=AirDensity;
%
% State space models can be generated at a single speed or multiple speeds.
% For flutter analysis we need a fine enough array of speeds.
%
% For flutter analysis (by another code) nspeeds speeds are used
% to span a range of speeds from SpeedLow to SpeedHigh.
% At each speed a state space model is
% generated, to be used later for flutter analysis by finding the
% eigenvalues of the [Ap] matrix as functions of speed.
%
nspeeds=10; %                                    <--------------------------nspeeds
%
%
% AIR SPEED interval over which to compute the state-space matrices
%
speedLow = 3; % m/s
speedHigh = 30; % m/s
%
% DYNAMIC PRESSURES corresponding to the above speeds
%
DynPressLow = 0.5*rho*speedLow*speedLow; % Pa
DynPressHigh = 0.5*rho*speedHigh*speedHigh; % Pa

if(nspeeds > 1)
    %dspeed=(speedHigh-speedLow)/(nspeeds-1);
    dDynPressure=(DynPressHigh-DynPressLow)/(nspeeds-1);
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
nstatesTotal=2*ns+ns*nLag+nLagG;
nstates=nstatesTotal;
%
% an array for storing the speeds used
%
DynamicPressures=linspace(DynPressLow,DynPressHigh,nspeeds)';
Speeds=sqrt(DynamicPressures*2./rho);

%
% Eqs. (1.38)
%
% an array storing the [Ap] matrices for the different speeds
ApMatrices=zeros(nstates,nstates,nspeeds);
% an array storing the [Bc] matrices for the different speeds
BcMatrices=[];
if(nc > 0)
    BcMatrices=zeros(nstates,3*nc,nspeeds);
end
% an array storing the [BG] matrices for the different speeds
BGMatrices=[];
if(igust == "YesGusts")
    BGMatrices=zeros(nstates,3,nspeeds);
end
%
% Array for storing values of highest frequencies (rad/sec) for which the
% math models are valid
%
OmegaMax=zeros(nspeeds,1);
%
%
for ispeed=1:nspeeds                                %     <------ Loop on Dyn Pressures
    %
    % Uspeed is in units consistent with the units used for analysis.
    % SpeedUnitConversionFactor convert speeds to units desired for
    % presentation of results
    %
    DynPress=DynamicPressures(ispeed);
    %
    Uspeed=Speeds(ispeed);
    %
    % Dynamic Pressure
    %
    %qD(ispeed)=0.5*rho*Uspeed^2;
    qD(ispeed)=DynPress;
    %
    %
    % Find the highest frequncy (rad / sec) for which the model is valid:
    %
    OmegaMax(ispeed)=kMax*Uspeed/b;
    %
    %
    % create Mbarbarss,Cbarbarss, and Kbarbarss (equations 1.17, 1.18)
    % only the ns x ns part is needed
    %
    Mbarbarss=zeros(ns,ns);
    Cbarbarss=zeros(ns,ns);
    Kbarbarss=zeros(ns,ns);
    %
    Mbarbarss(:,:)=Mss(:,:)-0.5*rho*b*b*Pbarss(:,:,3);
    Kbarbarss(:,:)=Kss(:,:)-0.5*rho*Uspeed*Uspeed*Pbarss(:,:,1);
    Cbarbarss(:,:)=Css(:,:)-0.5*rho*Uspeed*b*Pbarss(:,:,2);
    %
    Mbarbarsc=[];
    Cbarbarsc=[];
    Kbarbarsc=[];
    %
    if (nc > 0)
        %
        % create Mbarbarsc,Cbarbarsc, and Kbarbarsc (equations 1.17, 1.18)
        % only the NS x NC part is needed
        %
        %
        Mbarbarsc=zeros(ns,nc);
        Cbarbarsc=zeros(ns,nc);
        Kbarbarsc=zeros(ns,nc);
        %
        Mbarbarsc(:,:)= Msc-0.5*rho*b*b*Pbarsc(:,:,3);
        Kbarbarsc(:,:)= -0.5*rho*Uspeed*Uspeed*Pbarsc(:,:,1);
        Cbarbarsc(:,:)= -0.5*rho*Uspeed*b*Pbarsc(:,:,2);
        % Note: no [Ksc] and [Csc] here.
        %
        %
        % end if loop on whether there are control surfaces (if nc ne 0)
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
    % Create the aerodynamic lag equations for the structural and control
    % motions (eq. 1.23, 1.24)
    %
    Ar=[];
    Brs=[];
    Brc=[];
    if(nLag > 0)
        %
        Ar=zeros(ns*nLag,ns*nLag);
        Brs=zeros(ns*nLag,ns);
        if(nc > 0);
            Brc=zeros(ns*nLag,nc);
        else
        end
        %
        for iLags=1:nLag;
            istart=(iLags-1)*ns+1;
            iend=iLags*ns;
            %
            for i=istart:iend;
                Ar(i,i)=-beta(iLags);
            end
            %
            for i=1:ns;
                for j=1:ns;
                    Brs(istart-1+i,j)=Pbarss(i,j,3+iLags);
                end
            end
            %
            if (nc > 0)
                for i=1:ns;
                    for j=1:nc;
                        Brc(istart-1+i,j)=Pbarsc(i,j,3+iLags);
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
        % Create the aerodynamic lag equations for the gust column (eqs. 1.25, 1.26)
        %
        if(nLagG > 0)
            ArG=zeros(nLagG,nLagG);
            BrG=zeros(nLagG,1);
            %
            for i=1:nLagG
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
            for i=1:ns;
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
            for i=1:ns;
                for j=1:nLagG;
                    PGr(i,j)=PbarGs(i,1,j+3);
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
        T2r=-0.5*rho*Uspeed*Uspeed*InvMbarbarss*Ir; % < ------------------------------
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
    if (nc > 0)
        %
        % T2c
        %
        Temp=zeros(ns,3*nc);
        %
        for i=1:ns;
            for j=1:nc;
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
    % Total number of states = nstates
    %
    Ap=zeros(nstates,nstates);
    Bc=[];
    if (nc >0)
        Bc=zeros(nstates,3*nc);
    else
    end
    BG=[];
    if(igust == "YesGusts")
        BG=zeros(nstates,3);
    else
    end
    %
    % The Ap matrix
    % -------------
    %
    for i=1:ns;
        Ap(i,ns+i)=1.;
    end
    %
    for i=1:ns;
        for j=1:ns;
            Ap(i+ns,j)=T21(i,j);
            Ap(i+ns,j+ns)=T22(i,j);
        end
    end
    %
    if( nLag > 0)
        for i=1:ns;
            for j=1:ns*nLag;
                Ap(i+ns,2*ns+j)=T2r(i,j);
            end
        end
    else
    end
    %
    %
    if (igust == "YesGusts")
        for i=1:ns;
            for j=1:nLagG;
                Ap(i+ns,2*ns+nLag*ns+j)=T2Gr(i,j);
            end
        end
    else
    end
    %
    if(nLag > 0)
        for i=1:nLag*ns; %
            for j=1:ns;
                Ap(2*ns+i,ns+j)=Brs(i,j);
            end
        end
        %
        for i=1:nLag*ns;
            for j=1:nLag*ns;
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
        Bc=zeros(nstates,3*nc);
        %
        for i=1:ns;
            for j=1:3*nc;
                Bc(ns+i,j)=T2c(i,j);
            end
        end
        %
        if(nLag >0)
            for i=1:nLag*ns;
                for j=1:nc;
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
        BG=zeros(nstates,3);
        %
        for i=1:ns;
            for j=1:3;
                BG(ns+i,j)=T2G(i,j);
            end
        end
        %
        if(nLagG > 0)
            for i=1:nLagG;
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
    % The ASE model - the state equations:
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
%
% For open-loop and closed-loop control simulations, the full state space
% models are needed plus the output equations.
%
% -----------------------------
% -     Output Equations      -
% -----------------------------
%
% Eigenvactors (mode shape vectors)from the NASTRAN model:
% PHI_T N x (NS+NC)
% Each PHI_i is a vector
% {x1 y1 z1 r1 p1 q1 x2 y2 z2 r2 p2 q2 ...}')
% N is the number of degrees of freedom in the NASTRAN finite element model
% Each node in the FE model has 6 degrees of freedom (x,y,z,p,q,r)
%
load ../ASEInputData/PHI_T.mat; % < ---------------------------------------------PHI_T
PHI=PHI_T;
%
% The matrix of ns structural modes used (each mode shape is a column)
%
PHIs=PHI(:,1:ns);
%
% The nc control modes used (rigid body control rotations)
% Each control surface rigid-rotation mode is given as a vector of motions
% of the full N degrees of freedom in the FE model.
%
PHIc=[];
if(nc >0 )
    PHIc=PHI(:,NS+1:NS+nc);
else
end
%
% Accelerometers (ACC):
% PHI_Z_ACC is the T3 displacement at the accelerometers' locations
% PHI_X_ACC is the T1 displacement at the accelerometers' locations
%
% Loading Grid ID and their locations
load ../ASEInputData/GRID_ID_XYZ.mat; % < ---------------------------------------GRID_ID_XYZ
%
%
% 4-29-2021 Accelerometer locations from John Berg's measurements
GRID_ID_ACC=[2223,2423,3008]; % <------------------------------------------GRID_ID_ACC
%
% Finding locations ACC grid points
N_ACC=length(GRID_ID_ACC);
%
naccelerometers=N_ACC;
% < ---------------------------------------------------------------- naccelerometers
for i=1:N_ACC
    ind_ACC(i)=find(GRID_ID_ACC(i)== GRID_ID_XYZ);
end
%
% The coordinates of the measurement points
%
for i=1:N_ACC
    location=ind_ACC(i);
    xlocation(i,1)=i;
    xlocation(i,2)=GRID_ID_XYZ(location,2);
    xlocation(i,3)=GRID_ID_XYZ(location,3);
    xlocation(i,4)=GRID_ID_XYZ(location,4);
end
%
%
%< -------------------------------------------------------------------sensors locations
% ' Locations of Measurement Points '
xlocation
xxxxx(:)=xlocation(:,2);
yyyyy(:)=xlocation(:,3);
%
%
%
%
%
%
%
sensorlocations=zeros(N_ACC,3);
for iii=1:N_ACC
    sensorlocations(iii,1)=iii
end
sensorlocations(:,2)=xxxxx(:);
sensorlocations(:,3)=yyyyy(:);
sensorlocations
%
% Finding Z-displacement for each accelerometer at each mode
for i=1:N_ACC
    PHI_Z_ACC(i,1:NS+NC)=PHI((ind_ACC(i)-1)*6+3,1:NS+NC);
end
% Finding X-displacement for each accelerometer at each mode (
% for i=1:N_ACC
    % PHI_X_ACC(i,1:NS+NC)=PHI((ind_ACC(i)-1)*6+1,1:NS+NC);
% end
%
%
PHIsACC(1:N_ACC,1:ns)=PHI_Z_ACC(:,1:ns);
% PHIsACC(N_ACC+1:2*N_ACC,1:ns)=PHI_X_ACC(:,1:ns);
%
% Eqs. 1.64 - 1.73
%
% nstates=nstatesTotal=2*ns+ns*nLag+nLagG;
Tdisp=zeros(ns,nstates);
Tvel=zeros(ns,nstates);
%
for ii=1:ns
    Tdisp(ii,ii)=1.;
    Tvel(ii,ii+ns)=1.;
end
%
% Convert acceleration outputs to g's
%
% gravitational acceleration in m/s^2
gravitation = 9.80665; % m/s^2 < ----------------------------------------------gravity
%                                                                          Note Units
%
% Provide accelerations in g's
PHIsACC=PHIsACC/gravitation;
%
% Loop over speeds again
%
for ispeed=1:nspeeds
    %
    % {Vector of accelerations}=[Cp]{xp}+[Dcp]{uc}+DGp]{ug}
    %
    Ap(:,:)=ApMatrices(:,:,ispeed);
    %
    CpACC=PHIsACC*Tvel*Ap;
    %
    DcpACC=[];
    Bc=[];
    if(nc >0 )
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
    %
    % Vector of root loads
    %
    % Loads at the root FM(1:6,i) all 6 components for a mode i.
    % Z force for all modes is FM(3,:)
    load ../ASEInputData/FM.mat; % <-------------------------------------------------FM
    %                                                                    ---------Root
    % {RootLoads}=[FM}{qs}=[FM][Tdisp]{xp}
    %
    FMs=FM(:,1:ns); % dimensions (6) x (ns)
    RootMx=FMs(4,:); % dimensions (1) x (ns); x-moment at root
    RootStrain=-(RootMx*(0.0015875/2)/1.04016233e-7)/(6.89e10); % yStrain = Mx*(t/2)/(I*E)
    CpRootLoads=RootStrain*Tdisp;
    %
    % Expand the output equations to include loads after accelerations
    % nLoads loads outputs follow naccelerometers acceleration outputs
    %< ---------------------------------------------------------------------- nLoads
    nLoads=1;
    CpACC(naccelerometers+1:naccelerometers+nLoads,1:nstates)= ...
        CpRootLoads(:,:);
    if(nc >0)
        DcpACC(naccelerometers+1:naccelerometers+nLoads,1:3*nc)=0.;
    end

    if(igust=="YesGusts")
        DGpACC(naccelerometers+1:naccelerometers+1+nLoads,1:3)=0.;
    end

    %
    % Expand the output equations to include pitch sensor after loads
    % 
    % <---------------------------------------------------------------------- pitch

    GRID_ID_PITCH = [2003]; % wing root node
    idx_pitch_node = find(GRID_ID_XYZ(:,1)==GRID_ID_PITCH);
    idx_pitch_dof = ((idx_pitch_node-1)*6)+5; % 5th DOF b/c want Y-rotation
    PHI_PITCH = PHI(idx_pitch_dof,1:ns);

    nPitch = 2; % pitch and pitch derivative
    CpPitch = zeros(2,nstates); % <TODO>,<INSERT>
    CpPitch = [PHI_PITCH*Tdisp; PHI_PITCH*Tvel];
    CpACC(naccelerometers+nLoads+1:naccelerometers+nLoads+nPitch,1:nstates) = CpPitch;
    if(nc >0)
        DcpACC(naccelerometers+nLoads+1:naccelerometers+nLoads+nPitch,1:3*nc)=0.;
    end

    if(igust=="YesGusts")
        DGpACC(naccelerometers+nLoads+1:naccelerometers+nLoads+nPitch,1:3)=0.;
    end

    %
    % Total number of outputs
    %
    % < -------------------------------------------------------------- nTotalOutputs
    nTotalOutputs=naccelerometers+nLoads+nPitch;
    %
    %
    % The ASE model (p for "plant"):
    % {xpdot]=[Ap]{xp}+[Bc]{uc}+[BG]{uG}
    %
    % Accelerations
    % (yAcc}=[Cp]*{xp} + [Dpc]{uc}+[DpG]{uG}
    %
    % Root Loads
    % {Loads}=[CpRootLoads]{xp}
    %
    % {uc} =   {qc}
    %          s{qc}
    %         s^2{qc}
    %
    % {uG} = wG
    %       s wG
    %      s^2 wG
    %
    BcNewOrder=[];
    DcpACCNewOrder=[];
    %
    if(nc > 0)
        %
        % switch the [Bc] and [Dpc] matrices from working with:
        % {uc} =   {qc}
        %          s{qc}
        %         s^2{qc}
        % to working with:
        % {uc} =[ {qc1 s*qc1 s^2*qc1} {qc2 sqc2 s^2 qc2} ....}'
        %
        % control input order transformastion:
        %
        Tcontrolorder=zeros(3*nc,3*nc);
        %
        for ii=1:nc;
            location1=3*(ii-1)+1;
            location2=location1+1;
            location3=location1+2;
            Tcontrolorder(ii,location1)=1.;
            Tcontrolorder(ii+nc,location2)=1.;
            Tcontrolorder(ii+2*nc,location3)=1.;
        end
        %
        % change order of control input in {uc}
        BcNewOrder=Bc*Tcontrolorder;
        clear Bc;
        Bc=BcNewOrder;
        %
        DcpACCNewOrder=DcpACC*Tcontrolorder;
        clear DcpAcc;
        DcpACC=DcpACCNewOrder;
        %
    else
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??? needed at all ???
    BG=[];
    %
    if(igust == "YesGusts")
        BG(:,:)=BGMatrices(:,:,ispeed);
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%???
    %
    BpcMatrices(:,:,ispeed)=Bc(:,:);
    CpMatrices(:,:,ispeed)=CpACC(:,:);
    DpcMatrices(:,:,ispeed)=DcpACC(:,:);
    DpGMatrices(:,:,ispeed)=DGpACC(:,:);
    %
    % end loop over speeds
    %
end
%
%
% 
%
Sref=[]; % < -----------------------------------------------------------Sref (m2)
span=[]; % < ----------------------------------------------------------------span (m)
mac=[]; % < ------------------------------------------------------------------mac (m)
%
nACCoutputs=naccelerometers;
nLoadsOutputs=nLoads;
%
nconditions=nspeeds;
%
save('DataModel_3to30mps_15Modes_4Lags.mat', ...
    'ns', 'nc', 'nstates','nACCoutputs','nLoadsOutputs','nconditions', ...
    'Speeds', 'DynamicPressures','OmegaMax', ...
    'ApMatrices', 'BpcMatrices', 'BGMatrices', 'CpMatrices', 'DpcMatrices',...
    'DpGMatrices','AirDensity','Sref','span','b');
