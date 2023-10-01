clear all
close all

% original database from: 2023-04-17
% updated datanase from: 2023-08-01
    % this one missing FM.mat; I just made a copy of FM_calc.mat named FM.mat

%%
% Inputs:
%
% NS - number of input structural modes
NS=15;
%
% NC - number of input rigid body control surface modes
NC=4;
%
% nLag - number of lag terms in the Roger approximation 
% of the structural/control modes
nLag=4;
%
% nLagG - number of lag terms in the Roger approximation of the gust vector
%
nLagG=32;
%
% b - reference semichord
b=0.4/2;
%
% Reduced Frequencies: kbar
load ASEInputData/k_bar.mat; kbar=k_bar;
%
% kMax  - the highest reduced frequency used for Roger RFA fitting
kMax=max(kbar);
%
% NK - the size of reduced frequencies
NK=length(kbar);

% Pbar0, Pbar1, Pbar2, Pbar3, ... stored in a single array Pbar
% of dimensions (NS+NC) x (NS+NC) x (2+nLag)
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
% of dimensions (NS+NC) x (1) x (2+nLagG)
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
%           and control modes. Dimension nLag x 1
for i=1:nLag
    betabar(i)=1.7*kMax*i^2/(nLag+1)^2;
end
%
% betabarG - vector of lag terms used in the Roger RFA fit for the gust
%            vector. Dimension nLagG x 1
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
% Eigenvactors: PHI_T N x (NS+NC) (each PHI_i is transformed into a vector {x1 y1 z1 r1 p1 q1 x2 y2 z2 r2 p2 q2 ...}')
load ASEInputData/PHI_T.mat; PHI=PHI_T;
%
% zeta - vector of structural damping ratios for the structural modes. 
%        Dimension : NS x 1
zeta=zeros(NS,1);
%

% Nresponses - number of responses to be studied
%
% for each response - a vector ModalContributions including the modal 
%                     contributions to the response.
%                     For example, if the response is a vertical
%                     displacement at a structural node, then the vector
%                     ModalContributions would contain the displacements of
%                     the point in each of the mode shapes

%% Additional Matrices

% ss, cc, and sc matrices
M_ss=M(1:NS,1:NS);
M_cc=M(NS+1:end,NS+1:end);
M_sc=M(1:NS,NS+1:end);

K_ss=K(1:NS,1:NS);
K_cc=K(NS+1:end,NS+1:end);
K_sc=K(1:NS,NS+1:end);

for i=1:nLag+3
	Pbar_ss(:,:,i)=Pbar(1:NS,1:NS,i);
	Pbar_cc(:,:,i)=Pbar(NS+1:end,NS+1:end,i);
	Pbar_sc(:,:,i)=Pbar(1:NS,NS+1:end,i);
end

for i=1:nLagG+3
	PbarG_s(:,:,i)=PbarG(1:NS,1,i);
	PbarG_c(:,:,i)=PbarG(NS+1:end,1,i);
end


%% Exact aerodynamic and gust matrices if needed

% Aerodynamic Matrices: Q
load ASEInputData/QHH_T.mat; Q=QHH_T;

% Gust Vectors: QG
load ASEInputData/QHG_T.mat; QG=QHG_T;

for i=1:length(kbar)
	Q_ss(:,:,i)=Q(1:NS,1:NS,i);
	Q_cc(:,:,i)=Q(NS+1:end,1:NS+1:end,i);
	Q_sc(:,:,i)=Q(1:NS,1:NS+1:end,i);

	QG_s(:,:,i)=QG(1:NS,1,i);
	QG_c(:,:,i)=QG(NS+1:end,1,i);
end

%% Accelerometers (ACC): PHI_Z_ACC is the T3 displacement at the accelerometers' locations
% Loading Grid ID and their locations
load ASEInputData/GRID_ID_XYZ.mat;

% ACC Grid ID, the order is based on the numbering in the presentation: Bay 3, Bay 4, Bay 5, Bay 6, Tip LE, Tip TE)
GRID_ID_ACC=[];

% Finding locations ACC grid points
N_ACC=length(GRID_ID_ACC);
for i=1:N_ACC
    ind_ACC(i)=find(GRID_ID_ACC(i)==GRID_ID_XYZ);
end
 
% Finding Z-displacement for each mode
for i=1:N_ACC
	PHI_Z_ACC(i,1:NS+NC)=PHI((ind_ACC(i)-1)*6+3,1:NS+NC);
end


%% Actuators (ACTU)
% The actuation for wing will be translated into T3 at the TE (since these are plate elements and R5=0). The actuation for vanes is directly translated into R5

% ACTU Grid ID, the Flaps are TE and the Vanes are actuation points on EA (middle)
GRID_ID_ACTU=[2616 2620 3602 5001];

% Finding locations ACTU grid points
N_ACTU=length(GRID_ID_ACTU);
for i=1:N_ACTU
    ind_ACTU(i)=find(GRID_ID_ACTU(i)==GRID_ID_XYZ);
end


% Finding Z-displacement for each mode (besides pitch)
for i=1:N_ACTU-1
	PHI_Z_ACTU(i,1:NS+NC-1)=PHI((ind_ACTU(i)-1)*6+3,1:NS+NC-1);
end

% Finding R5-rotation for each mode (pitch only)
for i=N_ACTU:N_ACTU
	PHI_R5_ACTU(i,1:1)=PHI((ind_ACTU(i)-1)*6+5,NS+NC:NS+NC);
end



% Forces at the root FM(1:6,i) all 6 components for a mode i. Z force for all modes is FM(3,:)
load ASEInputData/FM.mat;

% Calculation Process

% % Root Grid ID
% % GRID_ID_FZ=1003;
% GRID_ID_FZ=[1003 2001];
% % Finding locations Root grid points
% N_FZ=length(GRID_ID_FZ);
% for i=1:N_FZ
%     ind_FZ(i)=find(GRID_ID_FZ(i)==GRID_ID_XYZ);
% end
% 
% clear KKG
% load Database/KGG_FIXED_15_BCFREE;
% 
% F=KGG*PHI;
% 
% for j=1:NS
% 	FM(j,1:6)=0;
% 	for ii=1:6
% 		for i=1:N_FZ
% 			FM(j,ii)=FM(j,ii)+F((ind_FZ(i)-1)*6+ii,j);
% 		end
% 	end	
% end
% 
% FM=FM';
% 
% save ASEInputData/FM_calc.mat FM

