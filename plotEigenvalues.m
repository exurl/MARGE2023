% Title: Plot Eigenvalues
% Author: Anthony Su
% Date: 2023-05-02
    % modified 2023-08-11

function [fig1,fig2] = plotEigenvalues(lambda,u,q,reLim,imLim)
% INPUTS
    % lambda : (N) x (nSpeeds) matrix of eigenvalues
    % q      : (nSpeeds) x (1) dynamic pressures
    % reLim  : real axis limits (optional)
    % imLim  : imag axis limits (optional)
% OUTPUTS
    % fig1   : figure of eig components
    % fig2   : figure of eigs on S-plane

nSpeeds = length(q);
for idxSpeed = 1:nSpeeds
    speedNames(idxSpeed) = string(['q=',num2str(q(idxSpeed))]);
end

%% plot eig components
% plotting setup
fig1 = figure;
tiledlayout(2,1);

% plotting data
for idxSpeed = 1:nSpeeds
    % x-axis (either u or q)
    % xData = q(idxSpeed)*ones(size(lambda(:,idxSpeed)),"like",lambda(:,idxSpeed));
    xData = u(idxSpeed)*ones(size(lambda(:,idxSpeed)),"like",lambda(:,idxSpeed));

    % plot eigenvalue real components
    nexttile(1)
    plot(xData,real(lambda(:,idxSpeed)),'k.')
    hold on

    % plot eigenvalue imaginary component
    nexttile(2)
    mask = ~isnan(real(lambda(:,idxSpeed)));
    plot(xData(mask),imag(lambda(mask,idxSpeed)),'k.')
    hold on
end

% plotting formatting
nexttile(1)
ylabel('$\mathrm{Re}(\lambda)$','Interpreter','latex')
if(exist('reLim','var'))
    ylim(reLim)
end
grid on

nexttile(2)
ylabel('$\mathrm{Im}(\lambda)$','Interpreter','latex')
% xlabel('$q_D$ (Pa)','Interpreter','latex')
xlabel('$u$ (m/s)','Interpreter','latex')
if(exist('imLim','var'))
    ylim(imLim)
end
grid on

%% plot eig on S plane
fig2 = figure;
baseColor = [1,0,0];
for idxSpeed = 1:nSpeeds
    % specify plot color
    color = baseColor - [1,0,0]*idxSpeed/nSpeeds + [0,0,1]*idxSpeed/nSpeeds;

    % plot eigs
    plot(real(lambda(:,idxSpeed)),imag(lambda(:,idxSpeed)),'.',"MarkerEdgeColor",color)
    hold on
end

if(~(exist('reLim','var') && exist('imLim','var'))) % NAND
    axis equal
end
if(exist('reLim','var'))
    xlim(reLim)
end
if(exist('imLim','var'))
    ylim(imLim)
end

legend(speedNames)
grid on
xlabel('$\mathrm{Re}(\lambda)$','Interpreter','latex')
ylabel('$\mathrm{Im}(\lambda)$','Interpreter','latex')

end