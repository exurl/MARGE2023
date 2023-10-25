% Title: Anthony ASE LTI SS Plant Generation
% Author: Anthony Su
% Date: 2023-08-03

function returnObj = Anthony_ASE_SS_plant_generation(nastranInputDir,NS,NC,ns,nc,nLag,nLagG,gusts,u,rho,b,zeta,aeroCorrections,flapCorrections,kNew)
% INPUTS
    % nastranInputDir : name of folder containing NASTRAN matrices
    % NS              : # structural states in NASTRAN data
    % NC              : # rigid control states in NASTRAN data
    % ns              : # structural states in model output
    % nc              : # rigid control states in model output
    % nLag            : # structure+control aero lag states
    % nLagG           : # gust aero lag states
    % gusts           : boolean, whether gusts are modeled
    % u               : air speeds (flight conditions)
    % rho             : air density
    % b               : ref. semi-span used to generate Roger matrices
    % zeta            : uncoupled damping ratios of structural modes
    % aeroCorrections : multiplier on columns of aero matrix (ns)x(ns)x(3+nLag) or (ns)x(ns)
    % flapCorrections : multiplier on control surface aero forces (1)x(nc)
    % kNew            : manually set [K] diagonal elements
% OUTPUTS
    % A        : plant dynamics matrices
    % Bc       : control input dynamics matrices
    % BG       : gust input dynamics matrices
    % q        : dynamics pressures for each SS model
    % u        : speeds for each SS model
    % omegaMax : max valid freq. for each SS model
    % nStates  : # states

%% INPUT VALIDATION

nastranInputDir = char(nastranInputDir);
assert(numel(NS)==1)
assert(numel(NC)==1)
assert(numel(ns)==1 && ns<=NS)
assert(numel(nc)==1 && nc<=NC)
assert(numel(nLag)==1)
assert(numel(nLagG)==1)
assert(numel(gusts)==1 && class(gusts)=="logical")
assert(numel(u)>=0)
assert(numel(rho)==1)
assert(numel(b)==1)
assert(numel(zeta)==ns)
assert(all(size(aeroCorrections)==[ns,ns,3+nLag]) || all(size(aeroCorrections)==[ns,ns]))
assert(numel(flapCorrections)==nc)
assert(numel(kNew)==ns+nc)

if(gusts==false)
    nLagG = 0;
end

%% IMPORT DATA

% total number of states
nStates = (2*ns)+(nLag*ns)+(nLagG);

% total number of modes
N = NS+NC;

% import NASTRAN structure matrices
M = load("ASEInputData\MHH_T.mat").MHH_T;
K = load("ASEInputData\KHH_T.mat").KHH_T;

% import NASTRAN aerodynamic influence matrices
kBar = load([nastranInputDir,'k_bar.mat']).k_bar;
Pbar(:,:,1) = load([nastranInputDir,'A0.mat']).A0;
Pbar(:,:,2) = load([nastranInputDir,'A1.mat']).A1;
Pbar(:,:,3) = load([nastranInputDir,'A2.mat']).A2;
if(nLag>0)
    D = load([nastranInputDir,'D.mat']).D; % (N) x (nLag*N)
    E = load([nastranInputDir,'E.mat']).E; % (nLag*N) x (N)
    D = reshape(D,N,N,nLag); % (N) x (N) x (nLag)
    E = pagetranspose(reshape(E',N,N,nLag)); % (N) x (N) x (nLag)
    Pbar(:,:,4:3+nLag) = pagemtimes(D,E); % (N) x (N) x (nLag+3)
end

if(gusts)
    PbarG(:,1,1) = load([nastranInputDir,'AG0.mat']).AG0;
    PbarG(:,1,2) = load([nastranInputDir,'AG1.mat']).AG1;
    PbarG(:,1,3) = load([nastranInputDir,'AG2.mat']).AG2;
    if(nLagG>0)
        DG = load([nastranInputDir,'DG.mat']).DG; % (N) x (nLag*N)
        EG = load([nastranInputDir,'EG.mat']).EG; % (nLag*N) x (N)
        DG = reshape(DG,N,1,nLagG); % (N) x (N) x (nLag)
        EG = pagetranspose(reshape(EG',1,1,nLagG)); % (N) x (N) x (nLag)
        PbarG(:,:,4:3+nLagG) = pagemtimes(DG,EG); % (N) x (N) x (nLag+3)
    end
end

% lag terms used in aerodynamic influence matrices
kbar = load([nastranInputDir,'k_bar.mat']).k_bar;
for idxLag = 1:nLag
    betabar(idxLag)=1.7*max(kbar)*idxLag^2/(nLag+1)^2;
end
for idxLag = 1:nLagG
    betabarG(idxLag)=1.7*max(kbar)*idxLag^2/(nLagG+1)^2;
end

%% STRUCTURE MODEL DEFINITION

% manual adjustment of stiffness diagonal elements
K(kNew>0) = kNew(kNew>0);

% decompose mass & stiffness matrices into structure & control components
Mss = M(1:ns,1:ns);
Msc = M(1:ns,NS+1:NS+nc);

Kss = K(1:ns,1:ns);
Ksc = zeros(ns,nc); % should be zero; force it to be so

% compute uncoupled undamped natural frequencies (eigenvalues)
omega = sqrt(diag(Kss)./diag(Mss)); % (ns)x(ns), rad/s

% construct uncoupled damping matrix
Css = diag(2*zeta.*omega.*diag(Mss)); % (ns)x(ns)
Csc = zeros(ns,nc); % (ns)x(nc)

%% AERODYNAMIC MODEL DEFINITION

% apply viscosity corrections
Pbar(1:ns,1:ns,:) = Pbar(1:ns,1:ns,:).*aeroCorrections; % element-wise scaling
Pbar(1:ns,NS+1:NS+nc,:) = Pbar(1:ns,NS+1:NS+nc,:).*flapCorrections; % column-wise scaling

% decompose Roger matrices into structure & control components
Pbarss = Pbar(1:ns,1:ns,:);
Pbarsc = Pbar(1:ns,NS+1:NS+nc,:);

if(gusts)
    PbarGs = PbarG(1:ns,:,:);
end

%% FLIGHT CONDITIONS DEFINITION

% create vector of flight conditions
nSpeeds = length(u);
q = 0.5*rho*u.*u;

% max valid frequency
omegaMax = max(kbar)*u/b;

%% GENERATE {xDot} = [A]{x} + [Bc]{uc} + [BG]{uG}

% initialize state-space matrices
A = zeros(nStates,nStates,nSpeeds);
Bc = zeros(nStates,3*nc,nSpeeds);
if(gusts)
    BG = zeros(nStates,3,nSpeeds);
end

% loop over flight conditions
for idxSpeed = 1:nSpeeds

    % compute coupled aeroelastic system matrices (Eq. 1.17)
    Mbarbarss = Mss-0.5*rho*b*b*Pbarss(:,:,3);
    Cbarbarss = Css-0.5*rho*u(idxSpeed)*b*Pbarss(:,:,2);
    Kbarbarss = Kss-0.5*rho*u(idxSpeed)*u(idxSpeed)*Pbarss(:,:,1);

    Mbarbarsc = Msc-0.5*rho*b*b*Pbarsc(:,:,3);
    Cbarbarsc = Csc-0.5*rho*u(idxSpeed)*b*Pbarsc(:,:,2);
    Kbarbarsc = Ksc-0.5*rho*u(idxSpeed)*u(idxSpeed)*Pbarsc(:,:,1);
    
    % lag roots from Roger fit to Laplace domain (Eq. 1.12)
    if(nLag>0)
        beta = betabar*u(idxSpeed)/b;
    end
    if(nLagG>0)
        betaG = betabarG*u(idxSpeed)/b;
    end

    % structure+control lag state dynamics (Eq. 1.23, 1.24)
    if(nLag>0)
        for idxLag = 1:nLag
            ArBlocks{idxLag} = -beta(idxLag)*eye(ns);
        end
        Ar = blkdiag(ArBlocks{:});
        
        BrsBlocks = num2cell(Pbarss(:,:,4:end),[1 2]);
        Brs = vertcat(BrsBlocks{:});

        BrcBlocks = num2cell(Pbarsc(:,:,4:end),[1 2]);
        Brc = vertcat(BrcBlocks{:});
    end

    % gust lag state dynamics (Eq. 1.25, 1.26)
    if(nLagG>0)
        ArG = diag(-betaG);
        BrG = ones(nLagG,1);
    end

    % [Ir] (Eq. 1.28, 1.29)
    if(nLag>0)
        Ir = repmat(eye(ns),1,nLag);
    end

    % [PGr] (Eq. 1.28, 1.30)
    if(nLagG>0)
        PGrBlocks = num2cell(PbarGs(:,:,4:end),[1,2]);
        PGr = horzcat(PGrBlocks{:});
    end

    % pre-compute inverse of [Mbarbarss]
    invMbarbarss = inv(Mbarbarss);

    % [T21] [T22] [T2r] [T2Gr] [T2c] [T2G] (Eq. 1.32)
    T21 = -invMbarbarss*Kbarbarss;
    T22 = -invMbarbarss*Cbarbarss;
    if(nLag>0)
        T2r = -q(idxSpeed)*invMbarbarss*Ir;
    else
        T2r = [];
    end
    if(nLagG>0)
        T2Gr = 0.5*rho*u(idxSpeed)*invMbarbarss*PGr;
    else
        T2Gr = [];
    end
    T2c = -invMbarbarss*[Kbarbarsc,Cbarbarsc,Mbarbarsc];
    if(gusts)
        T2G = 0.5*rho*u(idxSpeed)*invMbarbarss*[PbarGs(:,:,1),PbarGs(:,:,2),PbarGs(:,:,3)];
    end


    % [Ap] (Eq. 1.34)
    A(1:ns,ns+1:2*ns,idxSpeed) = eye(ns); % row 1
    A(ns+1:2*ns,:,idxSpeed) = [T21,T22,T2r,T2Gr]; % row 2
    if(nLag>0)
        A(2*ns+1:2*ns+ns*nLag,ns+1:2*ns+ns*nLag,idxSpeed) = [Brs,Ar]; % row 3
    end
    if(nLagG>0)
        A(2*ns+ns*nLag+1:end,2*ns+ns*nLag+1:end,idxSpeed) = ArG; % row 4
    end

    % [Bc] (Eq. 1.34)
    Bc(ns+1:2*ns,:,idxSpeed) = T2c; % row 2
    if(nLag>0)
        Bc(2*ns+1:2*ns+ns*nLag,nc+1:2*nc,idxSpeed) = Brc; % row 3
    end

    % [BG] (Eq. 1.34)
    if(gusts)
        BG(ns+1:2*ns,:,idxSpeed) = T2G; % row 2
        if(nLagG>0)
            BG(2*ns+ns*nLag+1:end,2,idxSpeed) = BrG; % row 4
        end
    end
end

% return object
returnObj.A = A;
returnObj.Bc = Bc;
if(gusts)
    returnObj.BG = BG;
end
returnObj.q = q;
returnObj.u = u;
returnObj.omegaMax = omegaMax;
returnObj.nStates = nStates;

end