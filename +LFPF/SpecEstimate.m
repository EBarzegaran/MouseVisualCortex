function [AS_Wave] = SpecEstimate(Data,freqs)

%% prepare the data for analysis

% reverse the electrode order, so the electrodes go from lowest to deepest layers
% But since bipolar is done before reversing the electrodes, should I negate the Data?
y           = Data.Y; 


% Set the initial parameters
%freqs       = 1:80;
Nodes       = size(y,2);
ElecLabels  = arrayfun(@(x) ['L' num2str(x)],1:Nodes,'uni',false);
Times       = Data.Times;

Y           = permute(y,[1 4 2 3]);
Y           = reshape(Y,size(Y,1)*size(Y,2),size(Y,3),size(Y,4));

%% Wavelet: for parameter estimations, we pull all the conditions together

% Wavelet transform of the data
S_Wave      = xwt_cmorl_nv(permute(Y,[3 1 2]), Data.srate, freqs, 1, 5);
% normalize, to minimize 1/f noise
S_WaveN     = (S_Wave-mean(S_Wave(:,:,:,(Times<0)&(Times>-.3)),4)); % Normalize Spectrum
% keep the autospectrum
AS_Wave     = reshape(S_WaveN,Nodes*Nodes,size(S_WaveN,3),size(S_WaveN,4));
AS_Wave     = AS_Wave(1:Nodes+1:end,:,:); %diagonal spectrum
% open up some space in the memory
clear S_WaveN

%% STOK
% ff      = .98; % the range of ff to search
% p      = 15;  % the range of p to search, which is equal to : 40 to 10 Hz
% srate   = Data.srate;
% % estimate AR parameters using STOK
% KF          = dynet_SSM_STOK(Y,p,ff);
% % estimate the spectrum using AR coefficients
% KF          = dynet_parpsd(KF,srate,freqs',2);
% S_STOK      = abs(KF.SS.^2);
% % normalize, to minimize 1/f noise
% S_STOKN     = S_STOK - mean(S_STOK(:,:,:,(Times<0)& (Times>-.3)),4);
% % keep the autospectrum
% AS_STOK     = reshape(S_STOKN,Nodes*Nodes,size(S_STOKN,3),size(S_STOKN,4));
% AS_STOK     = AS_STOK(1:Nodes+1:end,:,:); %diagonal spectrum
% % open up some space in the memory

end
