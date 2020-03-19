function [PDC,  ID_roi, Time, Freq] = ExtractRoiPDC(StokAll,ROI1,ROI2)
% retruns the PDC data between the two ROIs, if you are interested in the
% intra area connectivity use a same name for ROI1 and ROI2
% the ROI names should be from this set: {'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam'}
% OUTPUT:
    % PDC:      is a L1xL2xTxFxS, where L1 and L2 is the number of layers in ROI1 and ROI2, T is the number of time points,
                %F is the frequency bins and S is the number of animals with that pair of ROIs 

%%
    IDs = fieldnames(StokAll);
%%    
    for S = 1:numel(IDs)
        if S==1
            Time = StokAll.(IDs{S}).Times;
            Freq = StokAll.(IDs{S}).f;
        end
        
        ROIsizes = contains(StokAll.(IDs{S}).ROIs,'VIS')*3+3;
        ROIindices = [0 cumsum(ROIsizes)];
        idx1 = find(strcmpi(StokAll.(IDs{S}).ROIs,ROI1));
        idx2 = find(strcmpi(StokAll.(IDs{S}).ROIs,ROI2));
        if ~isempty(idx1) && ~isempty(idx2)
            PDC{S} = StokAll.(IDs{S}).PDC(ROIindices(idx1)+1:ROIindices(idx1+1),ROIindices(idx2)+1:ROIindices(idx2+1),:,1:numel(Time));
            ID_roi{S} = IDs{S};
        else
             ID_roi{S} = 'NaN';
        end
    end
    
%     Time = StokAll.(IDs{S}).Times;
%     Freq = StokAll.(IDs{S}).f;
    
    if exist('PDC','var')
        PDC = cat(5,PDC{:});
    else
        PDC = zeros(6,6,numel(Freq),numel(Time),0);
    end
    
end