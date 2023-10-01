% Title: Half-Power Damping
% Author: Anthony Su
% Date: 2023-07-11

function [zeta,f_max,f_1,f_2] = halfPowerDamping(f,FRF)
% DESCRIPTION
    % computes damping ratio from frequency response data using half-power 
    % bandwidth method. input must contain only the data around the 
    % response peak of interest.
% INPUTS
    % f     : frequency domain vector
    % FRF   : frequency response function vector
% OUTPUTS
    % zeta  : damping ratio
    % f_max : frequency at peak
    % f_1   : lower half-power frequency
    % f_2   : higher half-power frequency

% input validation
assert(length(f)==length(FRF),"f not same size as FRF")

% find max and half-power indices
[FRF_max,idx_max] = max(FRF);
FRF_half = FRF_max/sqrt(2); % half-power
idxs = find(diff(sign(FRF-FRF_half)));

% verify there's only two half-power frequencies
assert(length(idxs)<=2,"more than two half-power points found in data")
assert(length(idxs)>=2,"less than two half-power points found in data")
assert(sum(idxs>idx_max)==1,"peak not within half-power points")

% find max and half-power frequencies
f_max = f(idx_max);
f_1 = min(f(idxs));
f_2 = max(f(idxs));

% find damping ratio
zeta = (f_2-f_1)/(2*f_max);
end