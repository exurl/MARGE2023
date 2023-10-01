function[FRF,f,coh,Gxx_hat,Gyy_hat,Gxy_hat,Gyx_hat] = czt_FRF(isig,osig,N,w,x,fs)
%This function returns the frequency response function (FRF),
%the vector of frequencies, and the coherence using the chirp z transform.

%INPUTS:
%isig - input signal
%osig - output signal
%N - Window length.
%w = [f1 f2] - the min and max frequency for the chirp z transform
%x = 'on' or 'off' - printout the chirp-z transform parameters.
%fs - sampling frequency


%OUTPUTS:
%FRF - complex-valued vector corresponding to the DFT between f1 and f2
%f - frequencies for the DFT
%coh - real-valued vector corresponding to the coherence

%*******************************************************
%********** INPUT ERROR HANDLING ***********************
%*******************************************************

if length(isig) ~= length(osig)
    error('I/O signal lengths do not match. Cannot compute FRF.');
end

if nargin > 6
    error('Too many input arguments.')
end

%********************************************************
%input checking on N happens below.
%the default value for N is halfway between Nmin and Nmax
%********************************************************

if isempty(w)
    w = [0.1 3];%Hz -> default frequencies for czt  
end

%************************************************************
%range checking on the output flag, x, is accomplished below.
%************************************************************s

if isempty(fs)
    error('Sampling frequency not given. Cannot compute FRF.');
end

%*************************************************************
%********* INPUT RANGE CHECKING ******************************
%*************************************************************
if numel(isig) < 100
    error('Not enough I/O data. Cannot compute accurate FRF.')
end

% if numel(isig) > 5000000
%     error('I/O data vector too long. Cannot compute FRF.');
%     %note that 150000 corresponds to 5 minutes of data at 500 Hz.
% end

if numel(w) ~= 2
    error('Invalid w parameter. Cannot compute FRF.');
end

if w(1) > w(2)
    error('Lower input frequency is greater than upper input frequency. Cannot compute FRF.');
end

if w(2) > 300
    warning('Upper CZT frequency too high. Defaulting to 5 Hz.')
    w(2) = 5;
end

if w(1) < 0 || w(2) < 0
    error('Invalid sign on input frequencies. Cannot compute FRF.');
end

if strcmp(x,'on')||strcmp(x,'ON')
    output_flag = 1;
elseif strcmp(x,'off') || strcmp(x,'OFF')
    output_flag = 0;
elseif isempty(x)
    output_flag = 0;
else
    warning('Undefined value for printout parameters. Defaulting to no printout')
    output_flag = 0;
end

if fs <= 0
    error('Invalid sampling frequency. Cannot compute FRF.');
end
%**********************************************************
%************ END OF INPUT RANGE CHECKING *****************
%**********************************************************
%define total length of data
L = length(osig);

%define upper and lower frequencies
f1 = w(1);
f2 = w(2);

%find limits on the window length
%Nmin is created to make the window cover two periods of the lowest
%frequency. This value can become larger than Nmax if the lower/starting
%czt frequency is chosen too low compared to the length of isig/osig.
Nmin = floor(2*(1/f1)*fs);
Nmax = floor(L/2);

%==============================================
%==== SANITY CHECK ON f1 and Nmin ============
%==============================================
if Nmin > Nmax
    %this happens when the supplied starting czt-frequency requires too 
    %large of a window. Small frequency => large period => large window 
    %The solution is to raise the lower frequency and warn the user.
    f1 = round((8*fs/Nmax),2);%this ensures that Nmax > 4*Nmin
    warning('Starting CZT frequency too low for accurate resolution. Starting frequency raised.');
    fprintf('Starting CZT frequency changed to: %d\n',f1);    
    %recalculate Nmin
    Nmin = floor(2*(1/f1)*fs);
end
%==============================================

%input check on N
if isempty(N)
    N = floor((Nmax - Nmin)/2) + Nmin;
end

%===============================================
%=== for user-supplied window length N =========
%===============================================
if N > Nmax
    warning('Window length too great. Defaulting to Nmax.') 
    N = Nmax;
end

if N < Nmin
    warning('Window length too small. Defaulting to Nmin.');
    N = Nmin;
end
%================================================

%********************************************************

%define length of overlap for 80% of window length N
P = floor((8/10)*N);%80% overlap
%********* RANGE CHECK ON P *****************************
if P <= 1
    error('Cannot perform 80% overlap. Cannot compute FRF.');
end
%********************************************************

%=========================================================
%== Partition I/O data into matrix with columns as windows
%=========================================================
input = buffer(isig,N,P);
output = buffer(osig,N,P);
%no error checking on input or output at this time

%number of windows
nr = size(input,2);

%length of window
T = size(input,1);

%create the hann window
hann_window = hann(N);

%apply the hann window to each segment
input_hann = hann_window.*input;
output_hann = hann_window.*output;

%create the omega and starting parameters for the czt
omega = exp(-1i*2*pi*(f2 - f1)/(N*fs));
a = exp(1i*2*pi*f1/fs);

%perform the chirp z transform along each column
input_dft = czt(input_hann,N,omega,a);
output_dft = czt(output_hann,N,omega,a);

%create the frequency vector for the output
fn = (0:N-1)'/N;
f = (f2 - f1)*fn + f1;

%calculate the rough estimates
%Note that we're calculating Gxy instead of Gyx. Gxy will be used to
%calculate the H1 FRF, where the standard assumption is that the majority
%of noise is at the output.
Gxx_tilde = (2/T)*(abs(input_dft)).^2;
Gyy_tilde = (2/T)*(abs(output_dft)).^2;
Gxy_tilde = (2/T)*conj(input_dft).*output_dft;
Gyx_tilde = (2/T)*conj(output_dft).*input_dft;


%calculate the smooth spectral estimates
U = sqrt(3/8);
Gxx_hat = (1/(U*nr))*mean(Gxx_tilde,2);
Gyy_hat = (1/(U*nr))*mean(Gyy_tilde,2);
Gxy_hat = (1/(U*nr))*mean(Gxy_tilde,2);
Gyx_hat = (1/(U*nr))*mean(Gyx_tilde,2);




%calculate the H1 FRF
FRF = Gxy_hat./Gxx_hat; % this is complex valued.

%calculate the H2 FRF
% FRF = Gyy_hat./Gyx_hat;


%calculate the coherence
%coherence per Sergio is H1/H2
coh = (abs(Gxy_hat).^2)./(abs(Gxx_hat).*abs(Gyy_hat));%coherence

if output_flag == 1
    fprintf('========================================\n');
    fprintf('====== CHIRP Z TRANSFORM PARAMETERS ====\n');
    fprintf('Starting frequency: %d Hz\n',f1);
    fprintf('Ending frequency: %d Hz\n',f2);
    fprintf('Current window length: %d\n',N);
    fprintf('Minimum window length: %d\n',Nmin);
    fprintf('Maximum window length: %d\n',Nmax);
    fprintf('Number of windows: %d\n',nr);
    fprintf('========================================\n');
    fprintf('========================================\n');
end


