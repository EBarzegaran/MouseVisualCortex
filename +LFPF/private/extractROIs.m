function session = extractROIs(Data,probeinfo,LayerNames,LayerInfo,Summary,ROIs_all,session,savepath)
    

%% How to reduce LGN and LP

dopca = false;
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
        % extract layer indices
        [~,ind]=intersect(LayerNames,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
        Data.Y = TempData(:,LayerInfo(ind),:,:);
        session.(probeinfo.structure_acronyms{ROIs_ind(r)}) = Data;
        ROInames{p} = probeinfo.structure_acronyms{ROIs_ind(r)};
        p = p+1;
        %-----------------------Plot the results-----------------------
        FIG = figure;
        plot_lfp_expanded(Data.Y,Data.Times,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
        title('Original')

        set(FIG,'unit','inch','position',[5 5 8 4],'color','w')
        if ~exist(savepath)
            mkdir(savepath)
        end
        export_fig(FIG,fullfile(savepath,[probeinfo.structure_acronyms{ROIs_ind(r)}]),'-pdf');close
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
                
            else
                % extract indices
                [~,ind]     =   intersect(LayerNames,arrayfun(@(x) [probeinfo.structure_acronyms{ROIs_ind(r)} num2str(x)],1:3,'uni',false));
                % check if there is a nan index
                sumData    =   zeros(size(Y,1),3,size(Y,3),size(Y,4));
                Ind         =   LayerInfo(ind); Ind = Ind(end:-1:1);
                sumData(:,~isnan(Ind),:,:)    =   Y(:,Ind(~isnan(Ind)),:,:);
                Data.Y      = sumData;
                
                % 
            end
            
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
            export_fig(FIG,fullfile(savepath,[probeinfo.structure_acronyms{ROIs_ind(r)}]),'-pdf'); close
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




