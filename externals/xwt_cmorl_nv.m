
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function [S, Coh, iCoh]  = xwt_cmorl_nv(X, fs, freq, pad0, omega0)
%---------------------
% This function computes auto- & cross- spectra by using Morlet wavelet
% transform.
% 
% INPUT
% - X:          bivariate or multivariate signals in the form of 3D matrix [nTime, nTrial, nChannel]
% - fs:         data sampling rate in Hz
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - pad0:       0 or 1 (0 to do not pad with zeros, 1 to do padding with zeros)
% - omega0:     central frequency parameter in Morlet wavelet
% OUTPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs, nTimes]
% - Coh:        coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
% - iCoh:       imaginary-coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
%---------------------

[Nt, Ntr, Nc] = size(X);        % Nt: timepoints, Ntr: trials, Nc: channels 

S = zeros(Nc,Nc,length(freq),Nt);

for itrial = 1:Ntr
    for ichan = 1:Nc  %.........now taking the wavelet transform....
        Wx(:,:,ichan) = cwt_cmorl(X(:,itrial,ichan),fs,freq,pad0,omega0);
    end
    %........ computing spectral matrix ........
    for ii = 1:Nc
        for jj = 1:Nc
            s(ii,jj,:,:) = Wx(:,:,ii).*conj(Wx(:,:,jj));
        end
    end
    S = S + s;
    %WF(:,:,:,itrial) = Wx;
    %itrial
end
S = S/Ntr;      % dividing by the number of trials
clear X Wx s

%....... computing coherence .......
if nargout>1,
    for ii = 1:Nc
        for jj = 1:Nc
            Coh(ii,jj,:,:) = abs(S(ii,jj,:,:)) ./ sqrt(S(ii,ii,:,:).*S(jj,jj,:,:)) ;                % Coherence (Coh)
            if nargout>2,
                iCoh(ii,jj,:,:) = imag(  S(ii,jj,:,:) ./ sqrt(S(ii,ii,:,:).*S(jj,jj,:,:))  );       % imaginary-Coherence (iCoh)
            end;
        end
    end
    Coh  = permute(Coh,  [3 4 1 2]);
    if nargout>2, iCoh = permute(iCoh, [3 4 1 2]); end;
end;



%==========================================================================
function  wave = cwt_cmorl(signal, fs, freq, pad, omega0)
%---------------------
% This function computes the complex Morlet wavelet transform.
% (Used inside the function xwt_cmorl_nv.m)
% 
% INPUT
% - signal:     1-D timeseries
% - fs:         data sampling rate in Hz
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - pad:        0 or 1 (0 to do not pad with zeros, 1 to do padding with zeros)
% - omega0:     central frequency parameter in Morlet wavelet
% OUTPUT
% - wave:       complex Morlet wavelet transform
%---------------------
% Reference: Torrence & Compo, BAMS (1998).
%---------------------

n1      = length(signal);               % minimum scale: 2/fs
x(1:n1) = signal; %- mean(signal);

%-----zero-padding to speed up and to reduce edge-effects
if (pad == 1),
    x  = [x, zeros(1,1*n1)];
    n2 = length(x);
    x  = [x, zeros(1,2^nextpow2(n2)-n2)];   % nextpower of 2 zeros padded
end;
n  = length(x);
xf = fft(x);        % Fourier transform of the (padded) time series

%....construct wavenumber array used in transform 
k  = [1:fix(n/2)];
k  = k.*((2.*pi)/(n/fs));
k  = [0., k, -k(fix((n-1)/2):-1:1)];
k0 = omega0;                                        % central frequency parameter in Morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));      % Fourier-Factor (Morlet wavelet)

%....scale and wave arrays
scale  = 1./(fourier_factor*freq);
Nscale = length(scale);

wave = zeros(Nscale,n);     % define the wavelet array
wave = wave + 1i*wave;      % make it complex

% Loop through all scales and compute transforms
for iscale = 1:Nscale
    daughter = wavelet_kernel(k0,k,scale(iscale));	
    wave(iscale,:) = ifft(xf.*daughter);
end
wave = wave(:,1:n1);        % wavelet transform without the zero-padding




%==========================================================================
function daughter = wavelet_kernel(k0, k, scale)
%---------------------
% This function generates the daughter wavelets.
% (Used inside the function cwt_cmorl.m)
% 
% INPUT
% - k0:         central frequency parameter in Morlet wavelet
% - k:          vector of Fourier frequencies 
% - scale:      wavelet scale
% OUTPUT
% - daughter:   vector of the wave function
%---------------------
% Reference: Torrence & Compo, BAMS (1998), pp.65.
%---------------------

n = length(k);

expnt    = -(scale.*k - k0).^2/2.*(k > 0.);
norm     = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);

daughter = norm*exp(expnt);
daughter = daughter.*(k > 0.);                      % Heaviside step function




%==========================================================================
function S = extrapolate_for_zfreq(S)
%---------------------
% This function extrapolates 
% 
% INPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs, nTimes]
% OUTPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs+1, nTimes]
%               where the new frequncy vector is: freq1 = [0 freq]
%---------------------

[N1,~,~,Nt] = size(S);

for ii = 1:Nt
    for k1 = 1:N1
        for k2 = 1:N1
            Y  = log(squeeze(S(k1,k2,1:5,ii)));
            yi = interp1(1:5,Y,0:5,'linear','extrap');
            y(k1,k2,1,ii) = real(exp(yi(1)));
        end
    end
end
S = cat(3,y,S);




%==========================================================================
function causality = hz2causality(H, S, Z, fs)
%---------------------
% This function computes Granger-Geweke causality in the bivariate case.
% 
% INPUT
% - H:          transfer function
% - S:          matrix of auto- and cross-spectra [Nc, Nc, nFreqs]
% - Z:          noise covariance matrix
% - fs:         sampling frequency
% OUTPUT
% - causality:  Granger-Geweke causality
%---------------------

Nc = size(H,2);

for ii = 1:Nc
    for jj = 1:Nc
        if ii~=jj,
            zc    = Z(jj,jj) - Z(ii,jj)^2/Z(ii,ii);
            numer = abs(S(ii,ii,:));
            denom = abs(S(ii,ii,:)-zc*abs(H(ii,jj,:)).^2/fs);
            causality(jj,ii,:) = log(numer./denom);
        end;
    end
    causality(ii,ii,:) = 0;                 % self-causality set to zero
end
causality = permute(causality,[3 1 2]);     % [freq x channel from x channel to]


