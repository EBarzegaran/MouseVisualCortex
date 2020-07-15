function [S_STOK, freqs,Times,ROIs,KF,probeinfo] = STOKEstimate_All(session_Data,P,ff,Method,Orders,freqs,doPDC)


%%
if ~exist('doPDC','var')
    doPDC = true;
end
%% 1) Rearrange data
ROIs        = session_Data.ROIs;

%Orders = {'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam','VISmma'};

[C,ia,ib]  = intersect(Orders,ROIs);
[~,I] = sort(ia);
ROIs = C(I);
unitID = [];
ProbeID = [];
for roi = 1:numel(ROIs)
    y       = session_Data.(ROIs{roi}).Y;
    y       = permute(y,[1 4 2 3]);
    Y{roi}	= reshape(y,size(y,1)*size(y,2),size(y,3),size(y,4));
    Times   = session_Data.(ROIs{roi}).Times;
    srate   = session_Data.(ROIs{roi}).srate;
    unitID  = [unitID session_Data.(ROIs{roi}).unitID];
    ProbeID = [ProbeID repmat(str2double(session_Data.(ROIs{roi}).ProbeID),[1 numel(session_Data.(ROIs{roi}).unitID)])];
end

Y           = cat(2,Y{:});
% detrending and normalizing
% Y = Y./mean(mean(Y,1).^2,3);
% Y = Y - mean(mean(Y,1),3);
%Y = Y ./sqrt(sum(mean(Y,1).^2,3));
%% APPLY STOK
% final estimation of FC
KF          = dynet_SSM_STOK(squeeze(Y),P,ff);
if doPDC
    S_STOK      = dynet_ar2pdc(KF,srate,freqs,Method,1, 2,1);
else
    S_STOK = [];
end
%% probe data preparation
probeinfo.unitID = unitID;
probeinfo.probeID = ProbeID;
probeinfo.sessionID = str2double(session_Data.(ROIs{roi}).SessionID);
end

