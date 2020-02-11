function [PDCs, IDs, Time, Freq, ROIs] = ExtractAllRoiPDC(StokAll,ROIs)



if ~exist('ROIs','var') || isempty(ROIs)
    ROIs = {'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam'};
end

for roi1 = 1:numel(ROIs)
    for roi2 = 1:numel(ROIs)
        [PDCs{roi1,roi2},IDs{roi1,roi2},Time,Freq] = ExtractRoiPDC(StokAll,ROIs{roi1},ROIs{roi2});
    end
end

end