% Title: AeroServoElastic (ASE) State-Space Model Generation
% Author: Anthony Su
% Date: 3/15/2023

% Last edited 4/7/2023

% this is a refactor of "P2022_08_30ASEModelGeneration.m" by Eli Livne
% goal: generate state-space model s{x_p} = [A_p]{x_p} + [B_c]{u_c} + [B_G]{u_G}

% for LARGE wing

close all
clear
% clc
format long g

% reference geometry
b = 0.2; % ref. semi-chord; why not mac/2?? @Marat used this number for k
S_ref = 0.0838708;
span = 0.508;
mac = 0.1651;

%% Define Plant [M] and [K] Matrices
% no. of structural and control modes in NASTRAN data
NS = 15;
NC = 4;

% no. of structural and control modes used to generate model
ns = NS;
nc = NC;

% import generalized [M] and [K] matrices from NASRTAN
M = load("ASEInputData/MHH_t.mat").MHH_T; % (N)x(N)
K = load("ASEInputData/KHH_t.mat").KHH_T; % (N)x(N)

% decompose generalized [M] and [K] into structural and control components
M_ss = M(1:ns,1:ns);             % (ns) x (ns)
M_cc = M(NS+1:NS+nc,NS+1:NS+nc); % (nc) x (nc)
M_sc = M(1:ns,NS+1:NS+nc);       % (ns) x (nc)
M_cs = M_sc';                    % (nc) x (ns)
K_ss=K(1:ns,1:ns);               % (ns) x (ns)
    % note: K_sc, K_cs, K_cc assumed to be zero

% apply correction to K_ss(1) to force first bending mode freq = 2.3 Hz
K_ss(1,1) = 208.84083;

%% Define Plant Uncoupled Viscous [C] Matrix
% define viscous damping coefficients corresponding to structural modes
zeta = 0.015*ones(1,NS); % (1) x (NS)
    % note: structural/hysteretic damping g = 2*zeta

% uncoupled natural frequencies in rad/s
omega_uncoupled = sqrt(diag(K_ss)./diag(M_ss)); % (1) x (ns)

% uncoupled natural frequencies in Hz
freqs_uncoupled = omega_uncoupled/(2*pi); % (1) x (ns)
disp("Uncoupled Natural Frequencies")
disp(freqs_uncoupled)

% construct structural uncoupled [C] Matrix
C_ss = 2*diag(zeta(1:ns).*omega_uncoupled).*diag(diag(M_ss)); % (ns) x (ns)

%% Find Plant Structural Modes
% solve eigenvalue problem
lambda = eig(K_ss,M_ss); % (ns) x (1)

% natural frequencies in rad/s
omega = sqrt(lambda); % (ns) x (1)
omega = sort(omega);

% natural frequencies in Hz
freqs = omega/(2*pi); % (ns) x (1)
disp("Natural Frequencies (Hz):")
disp(freqs)

% highest natural frequency in rad/s
omega_highest = omega(end);
    % note: later, omega_max is used to denote max freq. for which the
    % model is valid. that is different from omega_highest.

%% Aerodynamics of Structure/Control Modes (Roger Matrices)
% number of lag terms in Roger approximation of structural/control modes
nLag = 4;

% import tabulated reduced frequencies
k_bar = load("ASEInputData\k_bar.mat").k_bar; % (1) x (NK)
NK = length(k_bar);

% maximum valid reduced frequency from data
k_max = max(k_bar);

% import structure+control Roger approximation matrices 1-3
P_bar(:,:,1) = load("ASEInputData\A0.mat").A0; % [P0]
P_bar(:,:,2) = load("ASEInputData\A1.mat").A1; % [P1]
P_bar(:,:,3) = load("ASEInputData\A2.mat").A2; % [P2]
    % note: P_bar has shape (N) x (N) x (nLag+3)

% import structure+control Roger approximation [P4], [P5], ... lag matrices
N = NS + NC;
if(nLag == 0)
    D = [];
    E = [];
elseif(nLag > 0)
    D = load("ASEInputData/D.mat").D; % (N) x (nLags*N)
    E = load("ASEInputData/E.mat").E; % (nLags*N) x (N)
    for i = 1:nLag
        D_ = D(:,(i-1)*N+1:i*N); % ith column block of [D], (N) x (N)
        E_ = E((i-1)*N+1:i*N,:); % ith row block of [E],    (N) x (N)
        P_bar(:,:,i+3) = D_*E_;  % [P_i] = [D_i]*[E_i]
    end
end

% extract structural component of structure/control Roger approximation matrices
P_bar_ss = P_bar(1:ns,1:ns,:); % (ns) x (ns)

% extract sc-component of structure/control Roger approximation matrices
if(nc == 0)
    P_bar_sc = [];
elseif(nc > 0)
    P_bar_sc = P_bar(1:ns,NS+1:NS+nc,:); % (ns) x (nc)
end

% aerodynamic correction factor due to thickness effect on cL_alpha
aeroCorrection = 0.92;
P_bar_ss = P_bar_ss*aeroCorrection;
P_bar_sc = P_bar_sc*aeroCorrection;

% aerodynamic corrections factors due to control surface aerodynamics
if(nc == 0)
    flapCorrections = [];
elseif(nc > 0)
    flapCorrections(1:3) = 0.7; % flap corrections
    flapCorrections(4) = 1.0; % (no correction) gust vane
    for i = 1:3
        P_bar_sc(:,:,i) = P_bar_sc(:,:,i)*diag(flapCorrections);
    end
end

% structure/control lag roots
if(nLag == 0)
    beta_bar = [];
elseif(nLag > 0)
    beta_bar = [1:nLag]*1.7*k_max/(nLag+1)^2; % (1) x (nLag)
    disp("beta bar:")
    disp(beta_bar')
end

%% Aerodynamics of Gusts (Roger Matrices)
% do we model gusts? (boolean)
haveGusts = false;

% number of lag terms in Roger appproximation of unsteady gusts
nLag_G = 32;

% override nLag_G if no gusts
if(haveGusts == false)
    nLag_G = 0;
end

% <insert> comment
if(haveGusts == false)
    P_bar_G = [];
elseif(haveGusts == true)
    % import gust Roger approximation matrices 1-3
    P_bar_G(:,:,1) = load("ASEInputData\AG0.mat").AG0; % [P0_G]
    P_bar_G(:,:,2) = load("ASEInputData\AG1.mat").AG1; % [P1_G]
    P_bar_G(:,:,3) = load("ASEInputData\AG2.mat").AG2; % [P2_G]
        % note: P_bar_G has shape (N) x (1) x (nLag_G+3)

    % import gust Roger approximation [P4_G], [P5_G], ... lag matrices
    if(nLag_G == 0)
        D_G = [];
        E_G = [];
    elseif(nLag_G > 0)
        D_G = load("ASEInputData/DG.mat").DG; % (N) x (nLag_g)
        E_G = load("ASEInputData/EG.mat").EG; % (nLag_g) x (1)
        for i = 1:nLag_G
            D_G_ = D_G(:,i); % (N) x (1)
            E_G_ = E_G(i,:); % (1) x (1)
            P_bar_G(:,:,i+3) = D_G_*E_G_; % P_i = {E_i}^T * {G_i}
        end
    end

    % extract structural effects component of gust Roger approximation lag matrices
    P_bar_G_s = P_bar_G(1:ns,1,:); % (ns) x (1) x (3+nLag_G)
    
    % aerodynamic correction factor due to thickness effect on cL_alpha
    P_bar_G_s = P_bar_G_s*aeroCorrection;

end

% gust lag roots
if(nLag_G == 0)
    beta_bar_G = [];
elseif(nLag_G > 0)
    beta_bar_G = [1:nLag_G]*1.7*k_max/(nLag_G+1)^2; % (1) x (nLag_G)
end

%% Flight Conditions
% air density in slug/in3
rho = 1.1464e-7;

% dynamic pressure range of interest in psf
q_D_low = 6;
q_D_high = 10;

% convert dynamic pressure range of interest to psi
q_D_low = q_D_low/12^2;
q_D_high = q_D_high/12^2;

% airspeed range of interest in in/s
nSpeeds = 3;
speedHigh = sqrt(2*q_D_high/rho);
speedLow = sqrt(2*q_D_low/rho);

%% Loop Over Flight Conditions
% linear range of dynamic pressure iterate over
if(nSpeeds > 1)
    q_D_vector = linspace(q_D_low,q_D_high,nSpeeds); % (1) x (nSpeeds)
end

% total number of states in model
nStates = (2*ns) + (ns*nLag) + (nLag_G);

% initialize data containers:
speeds = sqrt(q_D_vector*2/rho); % speeds
A_p_Matrices = zeros(nStates,nStates,nSpeeds); % [A_p]
if(nc==0)
    B_c_Matrices = [];
elseif(nc > 0)
    B_c_Matrices = zeros(nStates,3*nc,nSpeeds); % [B_c]
end
if(haveGusts == false)
    B_G_Matrices = [];
elseif(haveGusts == true)
    B_G_Matrices = zeros(nStates,3,nSpeeds); % [B_G]
end
omega_max = zeros(1,nSpeeds); % max frequencies @ which the model is valid
    % note: earlier, omega_highest is used to denote highest natural freq.
    % of the plant. that is different from omega_max.

% iterate over dynamic pressures
for iSpeed = 1:nSpeeds
    q_D = q_D_vector(iSpeed);

    % airspeed
    Uspeed = speeds(iSpeed);

    % max frequency for which the model is valid
    omega_max = k_max*Uspeed/b;

    % coupled aeroelastic [M_bar_ss], [K_bar_ss], [C_bar_ss]
    M_barbar_ss(:,:)= M_ss-0.5*rho*b*b*P_bar_ss(:,:,3);
    K_barbar_ss(:,:)= K_ss-0.5*rho*Uspeed*Uspeed*P_bar_ss(:,:,1);
    C_barbar_ss(:,:)= C_ss-0.5*rho*Uspeed*b*P_bar_ss(:,:,2);
        % note: plus [C_sc], but [C_sc] is [0]

    % coupled aeroelastic [M_bar_sc], [K_bar_sc], [C_bar_sc]
    if(nc == 0)
        M_barbar_sc = [];
        K_barbar_sc = [];
        C_barbar_sc = [];
    elseif(nc > 0)
        M_barbar_sc(:,:)= M_sc-0.5*rho*b*b*P_bar_sc(:,:,3);
        K_barbar_sc(:,:)= -0.5*rho*Uspeed*Uspeed*P_bar_sc(:,:,1);
            % note: plus [K_sc], but [K_sc] is [0]
        C_barbar_sc(:,:)= -0.5*rho*Uspeed*b*P_bar_sc(:,:,2);
            % note: plus [C_sc], but [C_sc] is [0]
    end

    % Convert beta_bar to beta
    if(nLag == 0)
        beta = [];
    elseif(nLag > 0)
        beta = Uspeed/b*beta_bar;
    end

    % convert beta_bar_G to beta_G
    if(nLag_G == 0)
        beta_G = [];
    elseif(nLag_G > 0)
        beta_G = Uspeed/b*beta_bar_G;
    end

    % structure/control lag state dynamics (eq. 1.23, 1.24)
    if(nLag == 0)
        A_r = [];
        B_rs = [];
        B_rc = [];
    elseif(nLag > 0)
        % initialize strcture/control matrices
        A_r = zeros(ns*nLag,ns*nLag);
        B_rs = zeros(ns*nLag,nc);
        if(nc == 0)
            B_rc = [];
        elseif(nc > 0)
            B_rc = zeros(ns*nLag,nc);
        end

        % construct structure/control matrices [A_r], [B_rs], [B_rc]
        for iLag = 1:nLag
            % A_r
            iRows = (iLag-1)*ns+1 : iLag*ns;
            iCols = (iLag-1)*ns+1 : iLag*ns;
            A_r(iRows,iCols) = eye(length(iRows))*-beta(iLag);

            % B_rs
            iRows = (iLag-1)*ns+1 : iLag*ns;
            iCols = 1:ns;
            B_rs(iRows,iCols) = P_bar_ss(:,:,3+iLag);

            % B_rc
            if(nc > 0)
                iRows = (iLag-1)*ns+1 : iLag*ns;
                iCols = 1:nc;
                B_rc(iRows,iCols) = P_bar_sc(:,:,3+iLag);
            end
        end
    end

    % gust lag state dynamics (eq. 1.25, 1.26)
    % construct gust matrices [A_rG], [B_rG]
    if(haveGusts == false)
        A_rG = [];
        B_rG = [];
    elseif(haveGusts == true)
        if(nLag_G == 0)
            A_rG = [];
            B_rG = [];
        elseif(nLag_G > 0)
            A_rG = eye(nLag_G)*-beta_G(i);
            B_rG = ones(nLag_G,1);
        end
    end

    % construct [I_r] (eq. 1.29)
    if(nLag == 0)
        I_r = [];
    elseif(nLag > 0)
        I_r = repmat(eye(ns),1,nLag);
    end
    
    % construct P_Gr (eq. 1.30)
    if(haveGusts == false)
        P_Gr = [];
    elseif(haveGusts == true)
        if(nLag_G == 0)
            P_Gr = [];
        elseif(nLag_G > 0)
            P_Gr = reshape(P_bar_G_s(:,:,4:end),[ns,nLag_G]);
        end
    end

    % plant STM components (eq. 1.32)
    inv_M_barbar_ss = inv(M_barbar_ss);

    T_21 = -inv_M_barbar_ss*K_barbar_ss;

    T_22 = -inv_M_barbar_ss*C_barbar_ss;

    if(nLag == 0)
        T_2r = [];
    elseif(nLag > 0)
        T_2r = q_D*inv_M_barbar_ss*I_r;
    end

    if(haveGusts == false)
        T_2Gr = [];
        T_2G = [];
    elseif(haveGusts == true)
        if(nLag == 0)
            T_2Gr = [];
        elseif(nLag > 0)
            T_2Gr = 0.5*rho*Uspeed*inv_M_barbar_ss*P_Gr;
        end
        T_2G = 0.5*rho*Uspeed*inv_M_barbar_ss*reshape(P_bar_G_s(:,:,1:3),[ns,3]);
    end

    if(nc == 0)
        T_2c = [];
    elseif(nc > 0)
        T_2c = -inv_M_barbar_ss*[K_barbar_sc,C_barbar_sc,M_barbar_sc];
    end

    % Construct plant STM from components (eq. 1.34)
    % [A_p]
    A_p = zeros(nStates,nStates);
    A_p(1:ns,ns+1:2*ns) = eye(ns); % block 1,2
    A_p(ns+1:2*ns,1:ns) = T_21; % block 2,1
    A_p(ns+1:2*ns,ns+1:2*ns) = T_22; % block 2,2
    A_p(ns+1:2*ns,2*ns+1:2*ns+ns*nLag) = T_2r; % block 2,3
    A_p(ns+1:2*ns,2*ns+ns*nLag+1:2*ns+nLag*ns+nLag_G) = T_2Gr; % block 2,4
    A_p(2*ns+1:2*ns+ns*nLag,ns+1:2*ns) = B_rs; % block 3,2
    A_p(2*ns+1:2*ns+ns*nLag,2*ns+1:2*ns+ns*nLag) = A_r; % block 3,3
    A_p(2*ns+1:2*ns+ns*nLag,2*ns+1:2*ns+ns*nLag) = A_rG; % block 3,3
    
    % [B_c]
    if(nc == 0)
        B_c = [];
    elseif(nc > 0)
        B_c = zeros(nStates,3*nc);
        B_c(ns+1:2*ns,:) = T_2c;
        if(nLag > 0)
            B_c(2*ns+1:2*ns+ns*nLag,nc+1:2*nc) = B_rc;
        end
    end

    % [B_G]
    if(haveGusts == false)
        B_G = [];
    elseif(haveGusts == true)
        B_G = zeros(nStates,3);
        B_G(ns+1:2*ns,:) = T_2G;
        if(nLag_G > 0)
            B_G(2*ns+ns*nLag+1:end,2) = B_rG;
        end
    end

    % store state-space matrices for this airspeed
    A_p_Matrices(:,:,iSpeed) = A_p;
    if(nc > 0)
        B_c_Matrices(:,:,iSpeed) = B_c;
    end
    if(haveGusts == true)
        B_G_Matrices(:,:,iSpeed) = B_G;
    end
end

clear A_p B_c B_G

%% Output Equations
% for "true accelerations" (not accounting for sensor dynamics)

% import mode shapes from NASTRAN
phi = load('ASEInputData/PHI_T.mat').PHI_T; % 

% extract control modes (rigid body rotations)
if(nc == 0)
    phi_c = [];
elseif(nc > 0)
    phi_c = phi(:,NS+1:NS+nc);
end

% import grid ID and locations from NASTRAN
grid_ID_xyz = load('ASEInputData\GRID_ID_XYZ.mat').GRID_ID_XYZ;
    % note: each column is a mode. every 6 entries in the column are
        % [dx, dy, dz, dxx, dyy, dzz] where dxx, dyy, dzz are rotations.
        % rotations are all zero in our simulation.

% accelerometer locations on grid
grid_ID_acc = [2223,2423];
    % note: correpsonding to John Berg measurement 4/29/2021
N_acc = length(grid_ID_acc);

% accelerometer indices in grid_ID_XYZ
idx_acc = zeros(1,N_acc);
for i = 1:N_acc
    idx_acc(i) = find(grid_ID_acc(i) == grid_ID_xyz);
end

% accelerometer locations in xyz space
location_acc = zeros(N_acc,4);
location_acc(:,1) = 1:N_acc; % first column is accelerometer #
location_acc(:,2:4) = grid_ID_xyz(idx_acc,2:4); % columns 2-4 are x,y,z

disp('Accelerometer Locations (#, X, Y, Z)')
disp(location_acc)


% for each mode get the z-displacement and x-displacement of accelerometers
for i = 1:N_acc
    phi_z_acc(i,1:N) = phi((idx_acc(i)-1)*6+3,1:NS+NC);
    phi_x_acc(i,1:N) = phi((idx_acc(i)-1)*6+1,1:NS+NC);
end

% extract structural (non-rigid body) component and concatenate
phi_z_acc_s = phi_z_acc(:,1:ns);
phi_x_acc_s = phi_x_acc(:,1:ns);
phi_acc_s = [phi_z_acc_s; phi_x_acc_s]; % z stacked on top of x

% [T_disp] (eq. 1.67)
T_disp = zeros(ns,nStates);
T_disp(1:ns,1:ns) = eye(ns);

% [T_vel] (eq. 1.69)
T_vel = zeros(ns,nStates);
T_vel(1:ns,ns+1:2*ns) = eye(ns);

% gravitational acceleration in in/s^2
gravity = 32.17405*12;

% accelerations in g's (?)
phi_acc_s = phi_acc_s/gravity;

% import root loads (6DOF) from NASTRAN
nLoads = 6;
FM = load('ASEInputData\FM.mat').FM; % (6) x (ns)
FM_s = FM(:,1:ns);

% [C_p_loads] (eq. 1.92)
C_p_loads = FM_s*T_disp;

for iSpeed = 1:nSpeeds
    % [A_p]
    A_p = A_p_Matrices(:,:,iSpeed);

    % [C] (eq. 1.76)
    C_p_acc = phi_acc_s*T_vel*A_p;

    if(nc==0)
        B_c = [];
        D_c_p_acc = [];
    elseif(nc > 0)
        % [B_c]
        B_c = B_c_Matrices(:,:,iSpeed);
        % [D_c]
        D_c_p_acc = phi_acc_s*T_vel*B_c;
    end

    if(haveGusts == false)
        B_G = [];
        D_G_p_acc = [];
    elseif(haveGusts == true)
        % [B_G]
        B_G = B_G_Matrices(:,:,iSpeed);
        % [D_G]
        D_G_p_acc = phi_acc_s*T_vel*B_G;
    end

    % expand [C_p_acc] to include loads after accelerations
    C_p_acc(2*N_acc+1:2*N_acc+nLoads,1:nStates) = C_p_loads;

    % expand [D_p_acc]  to include loads after accelerations
    if(nc > 0)
        D_c_p_acc(2*N_acc+1:2*N_acc+nLoads,1:3*nc) = 0;
    end

    % expand [D_G_p_acc] to include loads after accelerations
    if(haveGusts == true)
        D_G_p_acc(2*N_acc+1:2*N_acc+nLoads,1:3) = 0;
    end

    % total no. outputs
    nOutputs = 2*N_acc+nLoads;

    % re-order input vector definition such that
    % {u_c} = {{q_c}; s{q_c}; s^2{q_c}}
    % {u_G} = {{q_G}; s{q_G}; s^2{q_G}}
    
    % re-ordered [B_c], [D_c] via a transformation matrix
    if(nc > 0)
        T_c_order = zeros(3*nc);
        for i = 1:nc
            T_c_order(i,3*(i-1)+1) = 1;
            T_c_order(nc+i,3*(i-1)+2) = 1;
            T_c_order(2*nc+i,3*(i-1)+3) = 1;
        end
        % [B_c]
        B_c = B_c*T_c_order;
        % [D_c]
        D_c_p_acc = D_c_p_acc*T_c_order;
    end

    % state-space matrices for this airspeed
    B_p_c_Matrices(:,:,iSpeed) = B_c;
    C_p_Matrices(:,:,iSpeed) = C_p_acc;
    D_p_c_Matrices(:,:,iSpeed) = D_c_p_acc;
    if(haveGusts == false)
        D_p_G_Matrices = [];
    elseif(haveGusts == true)
        D_p_G_Matrices(:,:,iSpeed) = D_G_p_acc;
    end
end

% reference values for output
N_acc_output = 2*N_acc;
N_loads_output = nLoads;

% save('ASE-LTI-SS_6to10psf_15Modes_4Lags.mat', ...
%     'ns', 'nc','nStates','N_acc_output','N_loads_output','nSpeeds', ...
%     'speeds','q_D_vector','omega_max', ...
%     'A_p_Matrices','B_p_c_Matrices','B_G_Matrices','C_p_Matrices','D_p_c_Matrices','D_p_G_Matrices',...
%     'rho','S_ref','span','b');