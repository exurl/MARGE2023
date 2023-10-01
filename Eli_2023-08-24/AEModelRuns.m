%
% Model Runs with AE state space model
% Eli 8/23/2023
%
clear all ;
close all;
%
load('marge1_initial_models.mat', ...
     'ns', 'nc', 'nstates','Nconditions', ...
     'Speeds', 'DynamicPressures','OmegaMax', ...
     'ApMatrices', 'BpcMatrices', 'CpMatrices','DpcMatrices',...
     'AirDensity');
%
nstates
ns
nc
AirDensity
%
ispeed=3
%
Ap=ApMatrices(:,:,ispeed);
Bpc=BpcMatrices(:,:,ispeed);
Cp=CpMatrices(:,:,ispeed);
Dpc=DpcMatrices(:,:,ispeed);
%
' Ap'
size(Ap)
Ap
'Bpc '
Bpc
size(Bpc)
'Cp '
Cp
size(Cp)
'Dpc'
size(Dpc)
Dpc
%
' -------------- '
%
% input vector
%
u=zeros(3*nc,1);
% {u} here is made of {qc} then s{qc} and then s^2{qc}
%
OmegaLow=0.; % lowest frequency for FRF in rad/sec
OmegaHigh=4.*2*pi; % highest frequency for FRF rad/sec
nfreqs=100; % no. of frequencies for plotting FRFs
dOmega=(OmegaHigh-OmegaLow)/(nfreqs-1); % step size in omega rad/sec
Omega=-dOmega;
%
% start loop on frequencies (omegas)
%
for ifreq=1:nfreqs
%
Omega=Omega+dOmega;
j=sqrt(-1);
jw=j*Omega;
%
% input of control surface number (out of 4 control surfaces)
% 
CSnumber=1; %           input due to control surface number < --------------
%
% control surface unit rotations in the astate space model are in rads.
% Here we convert the actual amplitude of the control surface input (given
% in degrees) to radians
%
InputAmplitude=1. ; % 1 deg
%
u(CSnumber)=pi/180.*InputAmplitude;
% s*qc 
u(CSnumber+4)=jw*u(CSnumber);
% s^2*qc
u(CSnumber+8)=jw*jw*u(CSnumber);
%
% To see the I/O due to ONLY qc, comment out the lines for s*qc and s^2*qc
%
% [A]=[jw*I-Ap]
%
A=-Ap;
% Add jw along the diagonal
for ii=1:nstates
A(ii,ii)=A(ii,ii)+jw;
end
% inverst [jw*I-Ap]
Ainv=inv(A);
% Find the states {x}
x=Ainv*Bpc*u;
% The outouts: 1st: the rigid body pitch. Then accelerometers 1-3. No
% strain yet.
y=Cp*x+Dpc*u;
%
% The id of the output to study
OutputID=2 ; % In this case: accelerometer 1 in g's < -----------------
%
out=y(OutputID);
if(OutputID == 1)
out=out*180/pi; % the first output is the root rotation. It needs conversion to degrees
end
gain(ifreq)=abs(out);
phase(ifreq)=angle(out);
f(ifreq)=Omega/2/pi;
end
figure(1)
plot(f,gain)
grid on
title(['Output Gain ' num2str(OutputID) ' due to Input ' num2str(CSnumber)])
%
figure(2)
plot(f,phase*180./pi)
grid on
title(['Output phase (deg) ' num2str(OutputID) ' due to Input ' num2str(CSnumber)])