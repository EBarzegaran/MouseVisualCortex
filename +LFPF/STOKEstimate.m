function [S_STOK, freqs,Times] = STOKEstimate(Probe_Data,P,ff,Method,freqs)
Y = Probe_Data.Y;
Times = Probe_Data.Times;

% final estimation of FC
y = Y(:,:,:,:);% All conditions with contrast .8
y = permute(y,[1 4 2 3]);
Y = reshape(y,size(y,1)*size(y,2),size(y,3),size(y,4));
KF          = dynet_SSM_STOK(squeeze(Y),P,ff);
S_STOK      = dynet_ar2pdc(KF,Probe_Data.srate,freqs,Method,1, 2,1);

end