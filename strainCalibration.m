% Title: Strain Calibration
% Author: Anthony Su
% Date: 2023-09-04

close all
clear all
clc

%% EXPERIMENT

% load data
load('windTunnel/2023 summer/old/20230717/forceAndStrain.mat')
    % strain    : 4-element cell array, each containing a 5001-element vector
    % time      : 5001-element vector
    % loadConds :

% truncate data taken before settling in
cutoffIdxs = [4100,2700,3300,2300];

% find average values
for idx = 1:4
    avgStrain(idx) = mean(strain{idx}(cutoffIdxs(idx):end));
end

%% PREDICTION IN CUSTOMARY UNITS

% beam parameters
b = 3; % in
h = 0.0625; % in
I = b*h^3/12; % in^4
c = -h/2; % in
E = 10e6; % psi

% Marat's adjustment of I
I_meters = 2.541e-11; % m^4
E_Pa = 6.89e10; % Pa
I = I_meters/0.0254^4; % in^4
E = E_Pa*0.000145038;

% bending moment
P = loadConds*0.2248089431; % convert load from N to lb
M = loadConds*21; % lb-in

predictedStrain = -M*c/(E*I);

% convert from strain to microstrain
predictedStrain = predictedStrain*1e6;

% flip sign to match experimental convention
predictedStrain = -predictedStrain;

%% COMPARE EXPERIMENT AND PREDICTION

% strains
figure
plot(loadConds,avgStrain,'*')
hold on
plot(loadConds,predictedStrain,'o')

xlim([-0.8,0.8])
grid on
legend({'experiment','prediction'})
xlabel('load (N)')
ylabel('microstrain')

% ratios of strains
ratios = avgStrain./predictedStrain

figure
plot(loadConds,ratios,'.')

xlim([-0.8,0.8])
ylim([0,1])
grid on
xlabel('load (N)')
ylabel('ratio of experiment to prediction')

% ratios of strains plotted against predicted strain
figure
plot(predictedStrain,ratios,'.')

ylim([0,1])
grid on
xlabel('predicted microstrain')
ylabel('ratio of experiment to prediction')

% ratios of CORRECTED strains
figure
plot(loadConds,avgStrain./(predictedStrain*0.75),'.')

xlim([-0.8,0.8])
ylim([0,1.25])
grid on
xlabel('load (N)')
ylabel('ratio of experiment to corrected prediction')

%% QUADRATIC FIT OF STRAIN RATIO VS STRAIN

P = polyfit(predictedStrain,ratios,2)

% Intuition: the Euler-Bernoulli beam theory along with constant 0.75 
% correction captures first-order effects. A quadratic correction will help
% capture second-order effects.
% 
% Issue: this is a nonlinear effect, so it doesn't fit in our LTI model.
% Thus, we will stick with the constant 0.75 multiplier for now.