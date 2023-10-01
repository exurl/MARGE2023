% Title: Simulate MARGE Response
% Author: Anthony Su
% Date: 2023-08-01

clear all
format shortg

%% IMPORTING
path = './ASE_SS.mat';
aseModel = load(path);
A = aseModel.A;
Bc = aseModel.Bc;
C = aseModel.C;
Dc = aseModel.Dc;

% nStates
    % 15 flexible modes
    % 15 flexible mode derivatives
    % 15 flexible mode aero lags (1)
    % 15 flexible mode aero lags (2)
    % 15 flexible mode aero lags (3)
    % 15 flexible mode aero lags (4)
    % gust vane aero lag (1-32)
% nInputs
    % 4 control inputs [ail1, ail2, elev, gust]
% nOutputs
    % 3 accelerometer outputs [2223,2423,3008]
    % 1 root strain output
    % 1 pitch output
    % 1 pitch output derivative

omegaMax = aseModel.omegaMax;
rho = aseModel.rho;
u = aseModel.u';
q = 0.5*rho*u.*u;

nSpeeds = length(u)
nStates = size(A,1)
nInputs = size(Bc,2)
nOutputs = size(C,1)

% I/O names
inputNames = ["aileron 1","aileron 2","elevator","gust vane"];

outputNames = ["wing LE acc","wing TE acc","tail acc","root strain","pitch","pitch first derivative"];

%% COMPUTE STABILITY

% % compute open-loop stability (eigenvalues)
% for idxSpeed = 1:length(u)
%     lambda(:,idxSpeed) = eig(A(:,:,idxSpeed),"vector");
% end
% [~,idxLambdaMax] = max(real(lambda));
% for idxSpeed = 1:nSpeeds
%     lambdaMax(idxSpeed) = lambda(idxLambdaMax(idxSpeed),idxSpeed);
% end
% 
% % print open-loop stability analysis
% disp('open-loop most unstable eigenvalues:')
% disp('       q (Pa)    idxLambda   Re(lambda)   Im(lambda)')
% disp([q,idxLambdaMax',real(lambdaMax)',imag(lambdaMax)'])

%% PLOT EIGENVALUES

% % ignore lambda with Re(lambda)<<0
% lambdaClose = lambda;
% lambdaClose(real(lambda)<-1) = NaN;
% 
% % plot
% reLim = [-1,0.1];
% imLim = [-10,10];
% plotEigenvalues(lambdaClose,u,q,reLim,imLim);

%% COMPUTE TIME-DOMAIN RESPONSE
% scenario: aileron 1 dwell at 1.45 Hz with magnitude 1 degree

% % time span
% tSpan = [0,20];
% 
% % sinusoidal control input
% % inFreq = 1.45; % Hz
% % input = @(t) [sin(t*inFreq*(2*pi));zeros(3,length(t))]; % not called u b/c u is speeds
% 
% % zero control input
% input = @(t) zeros(4,length(t));
% 
% % zero initial state
% % x0 = zeros(size(A,1),1);
% 
% % first-wing-bent initial state
% x0 = zeros(size(A,1),1);
% x0(2) = 1; % mode shape is in the negative direction, so need to flip it back
% 
% % derivative function
% odeFunc = @(t,x) A(:,:,1)*x+Bc(:,:,1)*input(t);
% 
% % integrate
% [tOut,xOut] = ode89(odeFunc,tSpan,x0); % high-order integrator for multiscale dynamics
% xOut = xOut';
% 
% % control history
% uOut = input(tOut');
% 
% % output history
% yOut = C(:,:,1)*xOut+Dc(:,:,1)*uOut;

%% PLOT TIME-DOMAIN RESPONSE

% fig = figure;
% tl = tiledlayout(2,1);
% xlabel(tl,'Time (s)')
% 
% % plot input
% ax = nexttile(1);
% plot(tOut,uOut)
% legend(inputNames)
% title(ax,'Inputs')
% ylabel('Rotation (deg)')
% grid on
% 
% % plot output
% ax = nexttile(2);
% plot(tOut,yOut)
% legend(outputNames)
% title(ax,'Outputs')
% % ylabel('Acceleration (g)')
% grid on

%% COMPUTE FREQUENCY RESPONSE

% create state-space object
sys = ss(A,Bc,C,Dc);

% bode of all inputs
for idxSpeed = 1:nSpeeds
    sys_ = sys(:,:,idxSpeed);
    % omega = logspace(-0.5,1.5,201); % Hz
    omega = linspace(0.5,2); % Hz
    [mag(:,:,:,idxSpeed),phase(:,:,:,idxSpeed),] = bode(sys_,omega*2*pi);
        % indices of mag and phase are (idxOutput,idxInput,:,idxSpeed)
end

%% SAVE FRF DATA

% convert mag+phase into complex FRF
FRF = mag.*exp(phase*1j);

% reformat matrices
FRF = permute(FRF,[3,2,1,4]); % reorder indices to (:,idxOutput,idxInput,idxSpeed)

% duplicate frequency vector for each I/O combo:
freq = zeros(length(omega),nInputs,nOutputs);
for idxIn = 1:nInputs
    for idxOut = 1:nOutputs
        freq(:,idxIn,idxOut) = omega;
    end
end

for idxSpeed = 1:nSpeeds
    FRF_ = FRF(:,:,:,idxSpeed);

    % create struct
    obj.path = path;
    obj.title = "State Space Model";
    obj.u = [];
    obj.y = [];
    obj.t = [];
    obj.rate = [];
    obj.freq = freq;
    obj.coh = [];
    obj.H1_FRF = FRF_;
    obj.H2_FRF = FRF_;
    obj.Hv_FRF = FRF_;

    % store in array
    dataObjs(idxSpeed) = obj;
end

% save array of FRF data structs along with corresponding speeds
save('FRF_ASE_SS.mat','dataObjs','q')
disp('saved FRF_ASE_SS.mat');

%% PLOT FRFS

% % convert magnitude from linear to dB
% mag = mag2db(mag);
% 
% % custom color order
% lineStyles = {'-',':'};
% colors = {'k','r','g','b','c','m'};
% tab_ = combinations(lineStyles,colors);
% styleOrder = string([[tab_.colors{:}]',[tab_.lineStyles{:}]']);
% 
% % initialize figure
% fig = figure;
% fig.Position = [0,0,1000,800];
% tl = tiledlayout(nOutputs,nInputs);
% xlabel(tl,'Frequency (Hz)','FontSize',18);
% ylabel(tl,'Magnitude (dB)','FontSize',18);
% 
% % for each output,
% for idxOut = 1:nOutputs
%     % for each input,
%     for idxIn = 1:nInputs
%         % initialize tile
%         ax = nexttile();
%         hold on
%           % set(ax,'XScale','log')
%         grid on
% 
%         % title top row only
%         if(idxOut==1)
%             title(ax,outputNames(idxOut),'FontSize',12);
%         end
% 
%         % ylabel first column only
%         if(idxIn==1)
%             ylabel(ax,inputNames(idxIn),'FontSize',12)
%         end
% 
%         % remove x tick labels on axes that are not at the bottom
%         if(ax.Layout.Tile<10)
%             ax.XTickLabel = {};
%         end
% 
%         % for each speed,
%         for idxSpeed = 1:nSpeeds
%             plot(omega,squeeze(mag(idxOut,idxIn,:,idxSpeed)),styleOrder(idxSpeed));
%         end
%     end
% end
% 
% % generate speed legend
% for idxSpeed = 1:nSpeeds
%     labels{idxSpeed} = ['u=',char(num2str(u(idxSpeed),'%.1f'))];
% end
% leg = legend(labels,'Orientation','horizontal');
% leg.Layout.Tile = 'north';
% 
% % export figure
% print(strcat('./responsePlots/margeResponse'),'-dpng','-r300')