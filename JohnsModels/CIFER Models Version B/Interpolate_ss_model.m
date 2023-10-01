%Interpolate and plot state space elements versus dynamic pressure
clear all
close all
clc

% data folder
folder = 'CIFER models';
data_folder = strcat('C:\Users\CIFER-MATLAB\Documents\system_ID\marge1_state_space_models\',folder,'\');
x = dir(data_folder);
dyn_press_index = 1;
plot_results = 1;
q_new = 235;
%% pull out all matrices
for jj = 1:length(x)
    if length(x(jj).name) > 2
        if x(jj).name(9) == 'q'
            %this means we have a valid state space model
            d = load(strcat(data_folder,'\',x(jj).name));
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
[B,I] = sort(dyn_pressure_vector);
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
%% plot the results
if plot_results == 1
    %A matrix elements
    A21 = squeeze(A_matrices(2,1,:));
    A22 = squeeze(A_matrices(2,2,:));
    A23 = squeeze(A_matrices(2,3,:));
    A24 = squeeze(A_matrices(2,4,:));
    A41 = squeeze(A_matrices(4,1,:));
    A42 = squeeze(A_matrices(4,2,:));
    A43 = squeeze(A_matrices(4,3,:));
    A44 = squeeze(A_matrices(4,4,:));
    
    figure(1)
    subplot(2,4,1)
    scatter(dyn_pressure_vector,A21,'filled')
    ylabel('row 2','FontWeight','bold')
    title('column 1')
    grid on

    subplot(2,4,2)
    scatter(dyn_pressure_vector,A22,'filled')
    title('column 2')
    grid on

    subplot(2,4,3)
    scatter(dyn_pressure_vector,A23,'filled')
    title('column 3')
    grid on

    subplot(2,4,4)
    scatter(dyn_pressure_vector,A24,'filled')
    title('column 4')
    grid on

    %4th row
    subplot(2,4,5)
    scatter(dyn_pressure_vector,A41,'filled')
    ylabel('row 4','FontWeight','bold')
    title('column 1')
    grid on

    subplot(2,4,6)
    scatter(dyn_pressure_vector,A42,'filled')
    title('column 2')
    grid on

    subplot(2,4,7)
    scatter(dyn_pressure_vector,A43,'filled')
    title('column 3')
    grid on

    subplot(2,4,8)
    scatter(dyn_pressure_vector,A44,'filled')
    title('column 4')
    grid on

    sgtitle(strcat(folder,': A matrix elements'))

    %B matrix elements
    B21 = squeeze(B_matrices(2,1,:));
    B22 = squeeze(B_matrices(2,2,:));
    B23 = squeeze(B_matrices(2,3,:));
    B24 = squeeze(B_matrices(2,4,:));
    B41 = squeeze(B_matrices(4,1,:));
    B42 = squeeze(B_matrices(4,2,:));
    B43 = squeeze(B_matrices(4,3,:));
    B44 = squeeze(B_matrices(4,4,:));

    %plot B matrix elements
    figure(2)
    subplot(2,4,1)
    scatter(dyn_pressure_vector,B21,'filled')
    ylabel('row 2','FontWeight','bold')
    title('column 1')
    grid on
    
    subplot(2,4,2)
    scatter(dyn_pressure_vector,B22,'filled')
    title('column 2')
    grid on

    subplot(2,4,3)
    scatter(dyn_pressure_vector,B23,'filled')
    title('column 3')
    grid on

    subplot(2,4,4)
    scatter(dyn_pressure_vector,B24,'filled')
    title('column 4')
    grid on

    %4th row
    subplot(2,4,5)
    scatter(dyn_pressure_vector,B41,'filled')
    ylabel('row 4','FontWeight','bold')
    title('column 1')
    grid on

    subplot(2,4,6)
    scatter(dyn_pressure_vector,B42,'filled')
    title('column 2')
    grid on

    subplot(2,4,7)
    scatter(dyn_pressure_vector,B43,'filled')
    title('column 3')
    grid on

    subplot(2,4,8)
    scatter(dyn_pressure_vector,B44,'filled')
    title('column 4')
    grid on

    sgtitle(strcat(folder,': B matrix elements'))
   
    %C matrix elements
    C11 = squeeze(C_matrices(1,1,:));
    C13 = squeeze(C_matrices(1,3,:));
    C21 = squeeze(C_matrices(2,1,:));
    C23 = squeeze(C_matrices(2,3,:));
    C31 = squeeze(C_matrices(3,1,:));
    C32 = squeeze(C_matrices(3,2,:));
    C33 = squeeze(C_matrices(3,3,:));
    C34 = squeeze(C_matrices(3,4,:));

    %plot C matrix elements
    figure(3)
    subplot(3,4,1)
    scatter(dyn_pressure_vector,C11,'filled')
    ylabel('row 1','FontWeight','bold')
    title('column 1')
    grid on
    
    subplot(3,4,3)
    scatter(dyn_pressure_vector,C13,'filled')
    title('column 3')
    grid on
    
    %2nd row
    subplot(3,4,5)
    scatter(dyn_pressure_vector,C21,'filled')
    ylabel('row 2','FontWeight','bold')
    title('column 1')
    grid on

    subplot(3,4,7)
    scatter(dyn_pressure_vector,C23,'filled')
    title('column 3')
    grid on

    %3rd row
    subplot(3,4,9)
    scatter(dyn_pressure_vector,C31,'filled')
    ylabel('row 3','FontWeight','bold')
    title('column 1')
    grid on

    subplot(3,4,10)
    scatter(dyn_pressure_vector,C32,'filled')
    title('column 2')
    grid on

    subplot(3,4,11)
    scatter(dyn_pressure_vector,C33,'filled')
    title('column 3')
    grid on

    subplot(3,4,12)
    scatter(dyn_pressure_vector,C34,'filled')
    title('column 4')
    grid on

    sgtitle(strcat(folder,': C matrix elements'))

    %D matrix elements
    D31 = squeeze(D_matrices(3,1,:));
    D32 = squeeze(D_matrices(3,2,:));
    D33 = squeeze(D_matrices(3,3,:));
    D34 = squeeze(D_matrices(3,4,:));

    %plot the D matrix elements
    %4th row
    figure(4)
    subplot(1,4,1)
    scatter(dyn_pressure_vector,D31,'filled')
    ylabel('row 3','FontWeight','bold')
    title('column 1')
    grid on

    subplot(1,4,2)
    scatter(dyn_pressure_vector,D32,'filled')
    title('column 2','FontWeight','bold')
    grid on

    subplot(1,4,3)
    scatter(dyn_pressure_vector,D33,'filled')
    title('column 3','FontWeight','bold')
    grid on

    subplot(1,4,4)
    scatter(dyn_pressure_vector,D34,'filled')
    title('column 4','FontWeight','bold')
    grid on

    sgtitle(strcat(folder,': D matrix elements'))


        
end

% sys = ss(Anew,Bnew,Cnew,Dnew);
