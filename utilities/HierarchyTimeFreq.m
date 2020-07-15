
function [Hf, Hccf, BestPerm] = HierarchyTimeFreq(STOK_avg,Method)

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

tic
for f = 1:100
    %if mod(f,10)==0,disp(['Freq = ' num2str(f)]);end
    for t = 1:numel(Time)
        M = squeeze((HscoreMR(:,:,f,t)));
        %M = M./nansum(M(:));
        [Hf(f,t,:),Hccf(f,t),BestPerm(f,t,:)] = HierarchyExtract(M,Method);
        
    end 
end
toc
%[Hf,Hccf, BestPerm] = HierarchyExtractTimeFreq(HscoreMR,Method);

end
