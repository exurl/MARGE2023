function sys = interpolate_ss_fn(q_new,model_folder)
%Inputs
%q_new = dynamic pressure in pascals. Must be greater than 68 and less than
%343

%Model_folder = folder where the CIFER models are stored. This variable
%should end in a backslash

if q_new < 68 || q_new > 343
    error('Desired dynamic pressure out of range. Must be in [68,343] Pa.');
end


% data folder
% folder = 'CIFER models';
% data_folder = strcat('C:\Users\CIFER-MATLAB\Documents\system_ID\marge1_state_space_models\',folder,'\');
x = dir(model_folder);
dyn_press_index = 1;


%% pull out all matrices
for jj = 1:length(x)
    if length(x(jj).name) > 2
        if x(jj).name(9) == 'q'
            %this means we have a valid state space model
            d = load(strcat(model_folder,'\',x(jj).name));
            A_matrices(:,:,dyn_press_index) = d.sys.A;
            B_matrices(:,:,dyn_press_index) = d.sys.B;
            C_matrices(:,:,dyn_press_index) = d.sys.C;
            D_matrices(:,:,dyn_press_index) = d.sys.D;
    
            %pull out the dynamic pressure
            dyn_pressure_vector(dyn_press_index) = str2double(x(jj).name(10:length(x(jj).name)-4));
            dyn_press_index = dyn_press_index + 1;
        end
    end
end

%% put all the data in ascending order
[~,I] = sort(dyn_pressure_vector);
dyn_pressure_vector = dyn_pressure_vector(I);
A_matrices = A_matrices(:,:,I);
B_matrices = B_matrices(:,:,I);
C_matrices = C_matrices(:,:,I);
D_matrices = D_matrices(:,:,I);

%% interpolate
%Using the plotting results, we can see that all parameters vary linearly
%(a few are not quite linear but it's probably close enough)
Anew = zeros(size(d.sys.A));
Bnew = zeros(size(d.sys.B));
Cnew = zeros(size(d.sys.C));
Dnew = zeros(size(d.sys.D));


%A matrix interpolation
for ii = 1:size(Anew,1)
    for jj = 1:size(Anew,2)
        Anew(ii,jj) = interp1(dyn_pressure_vector,squeeze(A_matrices(ii,jj,:)),q_new);
    end
end
%B matrix interpolation
for ii = 1:size(Bnew,1)
    for jj = 1:size(Bnew,2)
        Bnew(ii,jj) = interp1(dyn_pressure_vector,squeeze(B_matrices(ii,jj,:)),q_new);
    end
end
%C matrix interpolation
for ii = 1:size(Cnew,1)
    for jj = 1:size(Cnew,2)
        Cnew(ii,jj) = interp1(dyn_pressure_vector,squeeze(C_matrices(ii,jj,:)),q_new);
    end
end
%D matrix interpolation
for ii = 1:size(Dnew,1)
    for jj = 1:size(Dnew,2)
        Dnew(ii,jj) = interp1(dyn_pressure_vector,squeeze(D_matrices(ii,jj,:)),q_new);
    end
end
sys_plant = ss(Anew,Bnew,Cnew,Dnew);

%% Add Actuator Dynamics
%Define actuator model
w_actuator = 6*2*pi;
zeta_actuator = 0.8;

A_servo = [0 1;-(w_actuator)^2 -2*zeta_actuator*w_actuator];
B_servo = [0;w_actuator^2];
C_servo = [1 0];
D_servo = 0;

%Gust Vane Dynamics
A_GV = 0;
B_GV = 0;
C_GV = 0;
D_GV = 1;

A_actuation = blkdiag(A_servo,A_servo,A_servo,A_GV);
B_actuation = blkdiag(B_servo,B_servo,B_servo,B_GV);
C_actuation = blkdiag(C_servo,C_servo,C_servo,C_GV);
D_actuation = blkdiag(D_servo,D_servo,D_servo,D_GV);

sys_actuation = ss(A_actuation,B_actuation,C_actuation,D_actuation);

%% Combine plant and actuation
sys = series(sys_actuation,sys_plant);

