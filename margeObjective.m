% Title: Anthony ASE LTI SS System Analysis/Optimization Objective Function
% Author: Anthony Su
% Date: 2023-08-07

% converted from Anthony_ASE_SS_driver.m into a function on 10-24-2023

function [residual,magError,phaseError] = margeObjective(x,mpWeight,wantSave)
% INPUTS
    % x        : vector of tuning parameters
        % x(1) : bending natural frequency (Hz)
        % x(2) : pitching mode damping ratio
        % x(3) : bending mode damping ratio
        % x(4) : ail1 effectiveness multiplier
        % x(5) : ail2 effectiveness multiplier
        % x(6) : elevator effectiveness multiplier
        % x(7) : gust vane effectiveness multiplier
        % x(8:end) : column-wise elements of multipliers to [P0], then
        %            [P1], then [P2], etc. e.g. if [P] is 2x2x3, then
        %            x(8:end) has P111,P211,P121,P221,P112,P212,...
    % mpWeight : ratio of magnitude error to phase error weight in residual
    % wantSave : boolean whether to save model. If blank, does not save
% OUTPUTS
    % residual   : scalar measure of model's FRF error from truth (experiment)
    % magError   : magntiude error of each y/u/q combo
    % phaseError : phase error of each y/u/q combo

    %% INDEPENDENT VARIABLES/PARAMETERS

    persistent nastranInputDir NS NC ns nc zeta gusts b nLag nLagG rho q g accIds strainId pitchId
    
    if(isempty(nastranInputDir))
        % NASTRAN data folder
        nastranInputDir = 'ASEInputData/';
        
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
        aeroCorrections = [1,1; 1,1];
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

        %% USE SMALL, SLOW MODEL

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
    end

    %% INPUT VALIDATION
    assert(numel(x)==7+(3+nLag)*ns*ns)
    if(~exist('wantSave','var'))
        wantSave = char([]);
    end
    
    %% EXTRACT MODEL TUNING DESIGN VARIABLES
    
    omega2 = x(1);
    zeta1 = x(2);
    zeta2 = x(3);
    flapCorrections = x(4:7);
    aeroCorrections = reshape(x(8:end),2,2,[]);
    
    % ==== NATURAL FREQUENCY REPLACEMENT ====
    omegan = zeros(1,ns+nc); % Hz
    % pitching frequency
    % omegan(1) = 
    % bending frequency
    omegan(2) = omega2;
    
    % ==== DAMPING RATIO SCALING ====
    % pitching damping
    zeta(1) = zeta1;
    % bending damping
    zeta(2) = zeta2;
    
    % ==== AERODYNAMICS ====
    % static aero corrections
    % aeroCorrections =
    % control surface aero corrections
    % flapCorrections =
    
    %% INTERMEDIATE VARIABLES
    
    % stiffness matrix replacement
        % NOTE: omegan = sqrt(k/m) --> k = m*omegan^2, where m=1 in FEM
    kNew = omegan.^2;
    
    % flight conditions
    persistent u nSpeeds
    if(isempty(u))
        u = sqrt(2*q/rho);
        nSpeeds = length(u);
    end
    
    %% PLANT MATRICES GENERATION
    plantObj = Anthony_ASE_SS_plant_generation(nastranInputDir,NS,NC,ns,nc,nLag,nLagG,gusts,u,rho,b,zeta,aeroCorrections,flapCorrections,kNew);
    Ap = plantObj.A;
    Bpc = plantObj.Bc;
    if(gusts)
        BG = plantObj.BG;
    else
        BG = [];
    end
    
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
    sysAct(:,:,1:length(q)) = sysAct; % copy for each speed
    
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
    
    %% MODEL TUNING CHANGES (POST-GENERATION)
    
    % flip gust input sign
    Bc(:,4,:) = -Bc(:,4,:);
    Dc(:,4,:) = -Dc(:,4,:);
    
    % flip strain output sign
    C(4,:,:) = -C(4,:,:);
    Dc(4,:,:) = -Dc(4,:,:);
    
    %% EXPORT MATRICES
    if(wantSave)
        filename = "ASE_SS.mat";
        if(gusts)
            save(filename,'A','Bc','BG','C','Dc','DG','omegaMax','q','NS','NC','ns','nc','nLag','nLagG','b','aeroCorrections','flapCorrections','rho','g','accIds','strainId','pitchId')
        else
            save(filename,'A','Bc','C','Dc','omegaMax','q','NS','NC','ns','nc','nLag','nLagG','b','aeroCorrections','flapCorrections','rho','g','accIds','strainId','pitchId')
        end
        disp(['model generated and saved at ',char(filename)])
    end
    
    %% CALL RESPONSE ANALYSIS FUCTIONS
    
    % compute FRFs
    sys = ss(A,Bc,C,Dc);
    omegaVec = 0.4:0.05:2; % Hz
    dataObjs = margeComputeFRF(sys,q,omegaVec);

    % compute error from experiment
    [magError,phaseError] = margeCompareFRF(dataObjs);

    % compute residual
    error = magError*sqrt(mpWeight) + phaseError/sqrt(mpWeight);
    residual = norm(error,'fro');

    % easy-read error return
    % error = reshape([sum(abs(error),1)],6,4,6);

end