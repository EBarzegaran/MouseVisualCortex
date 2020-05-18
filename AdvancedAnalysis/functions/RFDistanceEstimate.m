function StokALL = RFDistanceEstimate(StokALL,Probe_all)
% For now it is only the distance from the centroid is being calculated
% but we can also add RF overlap, etc...
%% 
IDs = fieldnames(StokALL);
    for id = 1:numel(IDs)
        Pinfo  = StokALL.(IDs{id}).Probeinfo;
        clear  RFCoords Pval;
        for x = 1:numel(Pinfo.unitID)
            [RFCoords(x,:) Pval(x)] = extract_coords(Pinfo.unitID(x),Pinfo.probeID(x),Pinfo.sessionID,Probe_all);
        end
        StokALL.(IDs{id}).ProbeRFDist = squareform(pdist(RFCoords));
        StokALL.(IDs{id}).ProbeRFPval = Pval;
    end
end

function [Cent Pval] = extract_coords(unitid,probeid,sessionid,Probe_all)
    ind = (Probe_all.Unit_ID==unitid) & (str2num(Probe_all.Probe_ID)==(probeid)) & (str2num(Probe_all.Session_ID)==(sessionid));
    [Pval cind] = min(Probe_all.RF_PValue(ind,:));
    Cent = Probe_all.RF_Centroid(ind,:,cind);
end