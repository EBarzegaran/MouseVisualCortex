function [R,MSE] = compare_spectrum(S1,S2)
% Compares two time-frequency spectra S1 and S2
% S1 and S2 are Nodes x frequency x time: should have the same sizes
% Should be updated later with a better and more solution
%%
% figure,
% subplot(1,2,1)
% plot(squeeze(mean(S1,2))')
% 
% subplot(1,2,2)
% plot(squeeze(mean(S2,2))')
%%

CCs = (arrayfun(@(x) corr2(squeeze(S1(x,:,:)),squeeze(S2(x,:,:))),1:size(S1,1)));
MSE = sqrt(sum(((S1(:)/max(S1(:)))-(S2(:)/max(S2(:)))).^2));
R = mean(CCs);

end