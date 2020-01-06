function session = extractROIs(Data,probeinfo,LayerInfo,summary,ROIs_all,session)
    


%% bipolar re-referencing should be done from top layer to bottom->
% should be multiplied by -1, here we only use .8 contrast 
Y = -1*Data.Y(:,:,:,Data.cnd_info.contrast==.8);
        
%% (1) first find the ROIs of interest
probeinfo.structure_acronyms(cellfun(@(C) any(isnan(C(:))), probeinfo.structure_acronyms))={'NOT'};
ROIs_ind = find(ismember(probeinfo.structure_acronyms,ROIs_all));


%% For each ROI:
p = 1;
ROInames = [];
for r = 1:numel(ROIs_ind)
    
    % (2) extract the data session
        % (a) for visual ROIs, indicate layers using laterinfo and summary
    if strfind(probeinfo.structure_acronyms{ROIs_ind(r)},'VIS') 
        TempData = Y(:,probeinfo.intervals(ROIs_ind(r)+1):-1:probeinfo.intervals(ROIs_ind(r))+1,:,:);
        Data.Y = TempData(:,LayerInfo,:,:);
        session.(probeinfo.structure_acronyms{ROIs_ind(r)}) = Data;
        ROInames{p} = probeinfo.structure_acronyms{ROIs_ind(r)};
        p = p+1;
    else
    
        % (b) for LGv, LGd and LP, use summary to indicate layers
        labels = probeinfo.intervals(ROIs_ind(r)+1):-1:probeinfo.intervals(ROIs_ind(r))+2;
        % if the layer width is more or equal to 6 electrodes->(6*40)
        if numel(labels)>=6
            % PCA and keep the first 6 PCs
            TempData = Y(:,labels,:,:);
            Ytemp = (permute(TempData,[2 3 1 4]));
            Dims = size(Ytemp);
            Ytemp = reshape(Ytemp,Dims(1),prod(Dims(2:end)));
            [coeff,score,latent] = pca(Ytemp');
            S = permute(reshape(score',Dims),[3 1 2 4]);
            Data.Y = S(:,1:6,:,:);
            session.(probeinfo.structure_acronyms{ROIs_ind(r)}) = Data;
            
            ROInames{p} = probeinfo.structure_acronyms{ROIs_ind(r)};
            p = p+1;
        end
        
    end   

end

%% Add ROIs to the list of rois in the session

if isfield(session,'ROIs')
    for i = 1:numel(ROInames)
        session.ROIs{end+1} = ROInames{i};
    end
else
    session.ROIs = ROInames;
end


end