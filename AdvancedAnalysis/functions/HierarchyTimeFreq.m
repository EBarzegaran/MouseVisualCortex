
function [Hf, Hccf] = HierarchyTimeFreq(STOK_avg)

%% Hierarchical organization based on siegle paper

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
NROIs   = numel(STOK_avg.ROIs);
% average the scores

for roi1 = 1:NROIs
   ind1 = (roi1-1)*6+1:roi1*6;
   for roi2 = 1:NROIs
       if roi2 ~=roi1
           ind2 = (roi2-1)*6+1:roi2*6;
           HscoreMR(roi2,roi1,:,:) = nanmean(nanmean(PDC(ind2,ind1,:,:),1),2);
       end
   end
end
%
%HscoreMR = PDC;

% for t = 50:numel(Time)
%     M = (nanmean(HscoreMR(:,:,:,t),3));
%     %M = M./nansum(M(:));
%     [H(t,:),Hcc(t)] = HierarchyExtract(M);
% end 

for f = 1:100
    for t = 50:numel(Time)
        M = squeeze((HscoreMR(:,:,f,t)));
        %M = M./nansum(M(:));
        [Hf(f,t,:),Hccf(f,t)] = HierarchyExtract(M);
    end 
end

end
