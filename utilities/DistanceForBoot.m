function DistanceBoots = DistanceForBoot(StokALL,BootIDs)
% BootIDs is a matrix of animal IDs for bootstraps with size nboots x
% bootsize
%% load ROIS
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
load ROInames;
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
NROIs = numel(ROIs);
RFsiglev = .05;
%% Keep in mind to exclude the insignificant RFs
for b = 1:size(BootIDs,1)
    
    % 1 - Select that subset of animals
    IDs     = {BootIDs{b,:}};
    %IDs     = setdiff(IDs,{'S821695405'});
    % 2 - Average over the distances of the animal subset
    Dists   = zeros(NROIs,NROIs,6,size(BootIDs,2));
    RFDists = zeros(NROIs,NROIs,6,size(BootIDs,2));
    
    for id = 1:numel(IDs)
        nrois       = numel(StokALL.(IDs{id}).ROIs);
        roiind      = find(ismember(ROIs,StokALL.(IDs{id}).ROIs));
        for roi1    = 1:nrois
            IND1    = (roi1-1)*6+1:roi1*6;
            for roi2 = roi1+1:nrois
                IND2    = (roi2-1)*6+1:roi2*6;
                Dists(roiind(roi1),roiind(roi2),:,id) = diag(StokALL.(IDs{id}).ProbeDist(IND1,IND2));
                rfdist = diag(StokALL.(IDs{id}).ProbeRFDist(IND1,IND2));
                rfdist(StokALL.(IDs{id}).ProbeRFPval(IND1)>RFsiglev | StokALL.(IDs{id}).ProbeRFPval(IND2)>RFsiglev) = NaN; % exclude the insignificant RFs
                 RFDists(roiind(roi1),roiind(roi2),:,id) = rfdist;
            end
        end
    end
    Dists_all(:,:,:,b)      = nanmean(Dists+permute(Dists,[2 1 3 4]),4);
    RFDists_all(:,:,:,b)    = nanmean(RFDists+permute(RFDists,[2 1 3 4]),4);
    RFDists_Nans(:,:,:,b)    = sum(isnan(RFDists+permute(RFDists,[2 1 3 4])),4);
   
end
 % 3 - Prepare output
DistanceBoots.Dists     =   Dists_all;
DistanceBoots.RFDists   =   RFDists_all;
DistanceBoots.RFDists_nans   =   RFDists_Nans;
end