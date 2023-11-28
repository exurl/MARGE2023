% Title: Model and Plots from {x}
% Author: Anthony Su
% Date: 2023-11-27

close all
clear all
clc

% load model
load('ns2_nLag0_UNTUNED.mat')
if(exist('x0','var'))
    x = x0;
end

% go back to main directory
cd ..

% generate ASE_SS.mat
margeObjective(x,mpWeight,true,ns,nLag);

% generate FRF_ASE_SS.mat
margeResponse;

% generate FRF plots in ./responsePlots/
margeFreqExperiment;