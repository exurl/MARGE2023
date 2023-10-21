% Title: System Compare with John's
% Author: Anthony Su
% Date: 2023-10-18

% returns my current ASE_SS.mat model along with John's interpolated model
% in the same standard MATALB state-space system object

function res = sysCompareJohn()
% OUTPUTS
    % res : (1) x (2) cell array containing
        % 1 : my state-space system from ASE_SS.mat (12 state 4 in 3 out)
        % 2 : John's state-space system from CIFER (4 state 4 in 3 out)

%% JOHN MODEL
% add johnSys directory to path
johnDir = 'JohnsModels/CIFER MOdels Version B/';
path(johnDir,path)

% import john's model
johnSpeeds = [68 94 169 209 283 343];
for i = 1:6
    speed = johnSpeeds(i);
    johnSys{i} = Interpolate_ss_fn(speed,'C:\\Users\Anthony\Documents\MATLAB\aeroelastic research\MARGE 2023\JohnsModels\CIFER Models Version B\');    
end

% re-order of john's model
stateOrder = [1,3,2,4];
inputOrder = [2,3,1,4];
outputOrder = [3,1,2];
for i = 1:6
    sys = johnSys{i};
    A = sys.A(stateOrder,stateOrder);
    B = sys.B(stateOrder,inputOrder);
    C = sys.C(outputOrder,stateOrder);
    D = sys.D(outputOrder,inputOrder);
    johnSys{i} = ss(A,B,C,D);
end

% extract into standard format
temp = cat(3,johnSys{1},johnSys{2},johnSys{3},johnSys{4},johnSys{5},johnSys{6});
clear johnSys
johnSys = temp;
clear temp

%% MY MODEL
% import my model
mySys = load('ASE_SS.mat');

% truncate my model
outputOrder = [1,4,5];
mySys.C = mySys.C(outputOrder,:,:);
mySys.Dc = mySys.Dc(outputOrder,:,:);

% extract into standard object
temp = ss(mySys.A,mySys.Bc,mySys.C,mySys.Dc);
clear mySys
mySys = temp;
clear temp

%% RETURN
res{1} = mySys;
res{2} = johnSys;

end