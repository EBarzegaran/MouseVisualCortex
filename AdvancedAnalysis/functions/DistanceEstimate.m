function StokALL = DistanceEstimate(StokALL,Probe_all)
    
    IDs = fieldnames(StokALL);
    for id = 1:numel(IDs)
        Pinfo  = StokALL.(IDs{id}).Probeinfo;
        clear Coords;
        for x = 1:numel(Pinfo.unitID)
            Coords(x,:) = extract_coords(Pinfo.unitID(x),Pinfo.probeID(x),Pinfo.sessionID,Probe_all);
        end
        StokALL.(IDs{id}).ProbeDist = squareform(pdist(Coords));
    end
end

function Coords = extract_coords(unitid,probeid,sessionid,Probe_all)
    ind = (Probe_all.Unit_ID==unitid) & (str2num(Probe_all.Probe_ID)==(probeid)) & (str2num(Probe_all.Session_ID)==(sessionid));
    Coords = [Probe_all.AP_CCF(ind) Probe_all.ML_CCF(ind) Probe_all.DV_CCF(ind)];
end