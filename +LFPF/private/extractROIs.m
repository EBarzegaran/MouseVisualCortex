function [session, probeinfo_organized] = extractROIs(Data,probeinfo,LayerNames,LayerInfo,Summary,ROIs_all,session,savepath,savename,StimParams,Sessions_ID,Probes_ID)
    

%% How to reduce LGN and LP

dopca = false;
%% bipolar re-referencing should be done from top layer to bottom->
% should be multiplied by -1, here we only use .8 contrast 
Pnames = fieldnames(StimParams);
CondSel = cellfun(@(x) ismember(Data.cnd_info.(x),StimParams.(x)),Pnames,'uni',false);
CondSel = prod(cat(3,CondSel{:}),3);
Y = -1*Data.Y(:,:,:,CondSel==1);
        
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
        % extract ROI, based on the structure index
        ROI_ind = probeinfo.intervals(ROIs_ind(r)+1):-1:probeinfo.intervals(ROIs_ind(r))+1;
        TempData = Y(:,ROI_ind,:,:);
        
        % extract layer indices
        [~,ind]=intersect(LayerNames,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
        Data.Y = TempData(:,LayerInfo(ind),:,:);
        
        Layer_ind = ROI_ind(LayerInfo(ind));
        Labels = cellfun(@(x) [probeinfo.structure_acronyms{ROIs_ind(r)} '_' x],LayerNames(ind),'uni',false);
        probeinfo_organized = prepare_probeinfo(probeinfo, Layer_ind, Labels, Sessions_ID, Probes_ID,Summary);
        
        % prepare session data
        Data.unitID     = probeinfo.Coords.id(Layer_ind);
        Data.ProbeID    = Probes_ID;
        Data.SessionID  = Sessions_ID;
        session.(probeinfo.structure_acronyms{ROIs_ind(r)}) = Data;
        ROInames{p}     = probeinfo.structure_acronyms{ROIs_ind(r)};
        
        p = p+1;
        
        
        %-----------------------Plot the results-----------------------
        FIG = figure;
        plot_lfp_expanded(Data.Y,Data.Times,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
        title('Original')

        set(FIG,'unit','inch','position',[5 5 8 4],'color','w')
        if ~exist(savepath)
            mkdir(savepath)
        end
        export_fig(FIG,fullfile(savepath,[probeinfo.structure_acronyms{ROIs_ind(r)}  savename]),'-pdf');close
        %--------------------------------------------------------------
    else
    
        % (b) for LGv, LGd and LP, use summary to indicate layers
        labels = probeinfo.intervals(ROIs_ind(r)+1):-1:probeinfo.intervals(ROIs_ind(r))+2;
        % if the layer width is more or equal to 6 electrodes->(6*40)
        if numel(labels)>=6
            TempData = Y(:,labels,:,:);
            if dopca
                % PCA and keep the first 6 PCs
                Ytemp = (permute(TempData,[2 3 1 4]));
                Dims = size(Ytemp);
                Ytemp = reshape(Ytemp,Dims(1),prod(Dims(2:end)));
                [coeff,score,latent] = pca(Ytemp');
                sumData = permute(reshape(score',Dims),[3 1 2 4]);
                Data.Y = sumData(:,1:6,:,:);
                probeinfo_organized = [];
                Layer_ind = [];
                
            else
                % extract indices
                [~,ind]     =   intersect(LayerNames,arrayfun(@(x) [probeinfo.structure_acronyms{ROIs_ind(r)} num2str(x)],1:3,'uni',false));
                % check if there is a nan index
                sumData    =   zeros(size(Y,1),3,size(Y,3),size(Y,4));
                Ind         =   LayerInfo(ind); Ind = Ind(end:-1:1);
                sumData(:,~isnan(Ind),:,:)    =   Y(:,Ind(~isnan(Ind)),:,:);
                Data.Y      = sumData;
                
                %
                Layer_ind = ROI_ind(LayerInfo(ind));
                Labels = cellfun(@(x) [probeinfo.structure_acronyms{ROIs_ind(r)} '_' x],LayerNames(ind),'uni',false);
                probeinfo_organized = prepare_probeinfo(probeinfo, Layer_ind, Labels, Sessions_ID, Probes_ID,Summary);

            end
            
                    
            % prepare session data
            Data.unitID     = probeinfo.Coords.id(Layer_ind);
            Data.ProbeID    = Probes_ID;
            Data.SessionID  = Sessions_ID;
            
            session.(probeinfo.structure_acronyms{ROIs_ind(r)}) = Data;
            ROInames{p} = probeinfo.structure_acronyms{ROIs_ind(r)};
            p = p+1;
            %-----------------------Plot the results-----------------------
            FIG = figure;
            subplot(1,2,1);
            plot_lfp_expanded(TempData,Data.Times,Summary(labels,:))
            title('Original')
            
            subplot(1,2,2);
            plot_lfp_expanded(sumData,Data.Times,[])
            title('PCs')
            
            set(FIG,'unit','inch','position',[5 5 15 5],'color','w')
            if ~exist(savepath)
                mkdir(savepath)
            end
            export_fig(FIG,fullfile(savepath,[probeinfo.structure_acronyms{ROIs_ind(r)} savename]),'-pdf'); close
            %--------------------------------------------------------------
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

if ~exist('probeinfo_organized','var')
    probeinfo_organized = [];
end
end


%%

function plot_lfp_expanded(Y, Time, labels)
   % first average
   YM = squeeze(mean(mean(Y,1),4));
   % find maximum
   M = max(abs(YM(:)));
   % plot 
   plot(Time,YM'+(0:M:(M*(size(YM,1)-1))),'linewidth',1.5);
   xlim([-.2 .3]);
   if exist('labels','var') && ~isempty(labels)
   	L = labels(:,2);L = cat(1,L{:});
   else
       L = [];
   end
   set(gca,'ytick',0:M:(M*(size(YM,1)-1)),'yticklabel',L);
   
   % adjust labels

end


function [probeinfo_organized] = prepare_probeinfo(probeinfo, Layer_ind, Labels, Sessions_ID, Probes_ID,Summary)
    % organize the probe information as a table for saving it later
    Unit_ID             = probeinfo.Coords.id(Layer_ind)';
    RF_Conds            = repmat(probeinfo.RF_mapping.Conds,[numel(Unit_ID),1]);
    Session_ID          = repmat(Sessions_ID,[numel(Unit_ID),1]);
    Probe_ID            = repmat(Probes_ID,[numel(Unit_ID),1]);
    AP_CCF              = probeinfo.Coords.AP_CCF(Layer_ind)';
    ML_CCF              = probeinfo.Coords.ML_CCF(Layer_ind)';
    DV_CCF              = probeinfo.Coords.DV_CCF(Layer_ind)';
    RF_Area             = probeinfo.RF_mapping.Area(:,Layer_ind)';
    RF_Centroid         = permute(probeinfo.RF_mapping.Centroid(:,:,Layer_ind),[3 1 2]);
    RF_PValue           = probeinfo.RF_mapping.PValue(:,Layer_ind)';
    RF_Map              = permute(probeinfo.RF_mapping.Maps(:,:,:,Layer_ind),[4 1 2 3]);
    Summary_Atlas       = cat(1,Summary{Layer_ind,2});
    
    probeinfo_organized = table(Unit_ID,Labels, Probe_ID, Session_ID, AP_CCF, ML_CCF, DV_CCF, RF_Area, RF_Centroid, RF_PValue, RF_Map,RF_Conds,Summary_Atlas);
end




