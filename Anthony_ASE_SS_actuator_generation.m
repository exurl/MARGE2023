% Title: Anthony ASE LTI SS Actuator Dynamics Generation
% Author: Anthony Su
% Date: 2023-08-18

function returnObj = Anthony_ASE_SS_actuator_generation()

%% SERVO TRANSFER FUCTION TO STATE-SPACE

% MKS HV6130 servo actuator transfer function:
% G(s) = y(s)/u(s) = k/(s^2 + a1*s + a0) where
k = 1461; % used to be 1394, adjusted to fix DC gain = 1
a1 = 62.2;
a0 = 1461;

[A,B,C,D] = tf2ss([0,0,k],[1,a1,a0]);

%% CONSTRUCT STATE-SPACE WITH OUTPUT [y yDot yDDot]'

% s*y(s) = [C]*s*{x} = [C]*([A]*{x} + [B]*u)
Cd = C*A;
Dd = C*B;

% s^2*y(s) = [C]*s^2*{x} = [C]*([A]*s*{x} + B*u) = [C]*([A]*([A]*{x} + [B]*u) + [B]*u)
% simplified: s^2*y(s) = [C]*[A]*[A]*{x} + ([C]*[A]*[B] + [C]*[B])*u
Cdd = C*A*A;
Ddd = C*A*B+C*B;

C = [C;Cd;Cdd];
D = [D;Dd;Ddd];

%% GUST VANE DYNAMICS
% gust vane dynamics are basically perfect but with a 0.034 s delay from
% the actuator driver hardware/software

[num,den] = pade(0.034,2);
[AG,BG,CG,DG] = tf2ss(num,den);

CGd = CG*AG;
DGd = CG*BG;

CGdd = CG*AG*AG;
DGdd = CG*AG*BG+CG*BG;

CG = [CG;CGd;CGdd];
DG = [DG;DGd;DGdd];

%% CONSTRUCT FOR ALL ACTUATORS
% assuming u = [u1 u1Dot u1DDot u2 ...]'

returnObj.A = blkdiag(A,A,A,AG);
returnObj.B = blkdiag(B,B,B,BG);
returnObj.C = blkdiag(C,C,C,CG);
returnObj.D = blkdiag(D,D,D,DG);
end