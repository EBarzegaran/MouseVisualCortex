
function [PostStim, PreStim] = FFT_estimate(Data,WinLength)

% Data is a LFP matrix
% WinLength is the length of FFT window in second

if ~exist('WinLength','var')
    WinLength =[];
end

y           = Data.Y; 
Times       = Data.Times;
Y           = permute(y,[1 4 2 3]);
Y           = reshape(Y,size(Y,1)*size(Y,2),size(Y,3),size(Y,4));
%% estimate ffts poststim
if ~isempty(WinLength)
    L = WinLength * Data.srate;
else
    L = numel(TimeInd);
end

TimeInd  = find(Times>0,1):numel(Times)-1;  % Time window
% L        = numel(TimeInd);                  % Signal length
if ~isempty(WinLength)
    L = WinLength * Data.srate;
else
    L = numel(TimeInd);
end
FS       = Data.srate;                      % sampling rate
Yfft     = fft(Y(:,:,TimeInd),L,3)./L;     % FFT & normalize
Yfft     = (Yfft(:,:,1:L/2+1));                  % final normalization
Yfft(:,:,2:end-1) = 2* Yfft(:,:,2:end-1);
F        = FS*(0:(L/2))/L;                  % FFT frequencies
%F        = F(1:FS/2+1);

FFT   = squeeze(mean(Yfft));
AMP   = squeeze(mean(abs(Yfft)));

PostStim.FFT = FFT;
PostStim.AMP = AMP;
PostStim.F = F;
%% estimate ffts pretstim
TimeInd  = 1:find(Times>=0,1)-1;  % Time window
%L        = numel(TimeInd);                  % Signal length
FS       = Data.srate;                      % sampling rate
Yfft     = fft(Y(:,:,TimeInd),L,3)./L;     % FFT & normalize
Yfft     = (Yfft(:,:,1:L/2+1));                  % final normalization
Yfft(:,:,2:end-1) = 2* Yfft(:,:,2:end-1);
F        = FS*(0:(L/2))/L;                  % FFT frequencies
%F        = F(1:FS/2+1);

FFT   = squeeze(mean(Yfft));
AMP   = squeeze(mean(abs(Yfft)));

PreStim.FFT = FFT;
PreStim.AMP = AMP;
PreStim.F = F;
end