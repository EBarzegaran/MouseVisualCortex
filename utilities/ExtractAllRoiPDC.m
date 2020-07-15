function [PDCs, IDs, Time, Freq, ROIs,RFDist,RFPval] = ExtractAllRoiPDC(StokAll,ROIs,Distanalysis)



if ~exist('ROIs','var') || isempty(ROIs)
    ROIs = {'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam'};
end

if ~exist('Distanalysis','var') || isempty(Distanalysis)
    Distanalysis = false;
end

for roi1 = 1:numel(ROIs)
    for roi2 = 1:numel(ROIs)
        [PDCs{roi1,roi2},IDs{roi1,roi2},Time,Freq, RFDist{roi1,roi2},RFPval{roi1,roi2}] = ExtractRoiPDC(StokAll,ROIs{roi1},ROIs{roi2},Distanalysis);
    end
end

end