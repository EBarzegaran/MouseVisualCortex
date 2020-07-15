function RF = RFDataplot(StokALL,Probe_all,P_th,Path)

IDs     = fieldnames(StokALL);

%% load the colors and the ROIs
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
Cobj = LFPF.RColors();
Colors = Cobj.MatrixColors(ROIs,'SubColors');
load ROInames.mat

%%
FIG1    = figure(1);% scatter plot of the centroids
set(FIG1,'color','w','unit','inch','position',[2 9 20 2.5])

FIG2    = figure(2);% area distribution
set(FIG2,'color','w','unit','inch','position',[2 6 20 2.5])
AreaInfo = [];

FIG3    = figure(3);% dist hist plots  
RFdist = cell(numel(ROIs));
set(FIG3,'color','w','unit','inch','position',[2 6 20 20])

for id = 1:numel(IDs)
    % first find the significant RFs
    Pinfo  = StokALL.(IDs{id}).Probeinfo;
    for x = 1:numel(Pinfo.unitID)
        [RFCoords(x,:) Pval(x) Area(x) label{x}] = extract_coords(Pinfo.unitID(x),Pinfo.probeID(x),Pinfo.sessionID,Probe_all);
    end
    
    
    labels = cellfun(@(x) strsplit(x,'_'),label,'uni',false);
    labels = cat(1,labels{:});
    ROIInd = arrayfun(@(x) find(ismember(ROIs,labels{x,1})),1:size(labels,1));
    LayerInd = arrayfun(@(x) str2num(labels{x,2}(end)),1:size(labels,1));

    % 
    figure(1);
    for roi = 1:numel(ROIs)
        subplot(1,numel(ROIs),roi);hold on;box on; grid on;
        Inds = find((ROIInd==roi) & (Pval<P_th));
        scatter(RFCoords(Inds,1),RFCoords(Inds,2),40,squeeze(Colors(roi,3,:))','filled')
        xlim([.5 9.5])
        ylim([.5 9.5])
        set(gca,'xtick',1:2:9,'xticklabel',-40:20:40,'ytick',1:2:9,'yticklabel',-40:20:40)
        if roi==1
            xlabel('X')
            ylabel('Y')
        end
        title(ROI_names.(ROIs{roi}))
    end
    
    
    % add areas
    AreaInfo_temp = [Area; ROIInd; LayerInd;Pval];
    AreaInfo = [AreaInfo AreaInfo_temp];
    
    
    % calculate distances...
    
    for roi1 = 1:numel(ROIs)
        Inds1 = find((ROIInd==roi1) & (Pval<P_th));
        CO1 = RFCoords(Inds1,:);
        for roi2 = 1:numel(ROIs)
            Inds2 = find((ROIInd==roi2) & (Pval<P_th));
            CO2 = RFCoords(Inds2,:);
            dist_temp = squareform(pdist([CO1;CO2]));
            if ~isempty(Inds1) & ~isempty(Inds2)
                if roi1==roi2
                    temp = dist_temp(tril(ones(size(dist_temp)))>0);
                else
                    temp = (dist_temp(1:numel(Inds1),numel(Inds1)+1: numel(Inds1)+numel(Inds2)));
                    temp = temp(:);
                end
                RFdist{roi1,roi2} = [RFdist{roi1,roi2}; temp];
            end
            
        end
        
    end
    
end
AreaInfo = AreaInfo';
%%
figure(2);

for roi = 1:numel(ROIs)
    subplot(1,numel(ROIs),roi);hold on;box on; grid on;
    Inds = find((AreaInfo(:,2)==roi) & (AreaInfo(:,4)<P_th));
    pdSix = fitdist(AreaInfo(Inds,1),'Kernel');
    x = 0:50:2500;
    ySix = pdf(pdSix,x);
    plot(x,ySix,'k-','LineWidth',2,'color',squeeze(Colors(roi,3,:)));
    xlim([0 2500])
    if roi==1
        ylabel('PDF')
        xlabel('Area')
    end
    title(ROI_names.(ROIs{roi}))
    set(gca,'ytick',[])
    YLIMs = ylim;
    text(2000, YLIMs(2)*.9,['n=' num2str(numel(Inds))])
end

%%

figure(3);
for roi1 = 1:numel(ROIs)
    for roi2 = 1:roi1
        subplot(numel(ROIs),numel(ROIs),(roi1-1)*numel(ROIs)+roi2);hold on;box on; grid on;
        set(gca,'ytick',[])
        %histogram(RFdist{roi1,roi2}*10,20)
        xlim([0 90])
         pdSix = fitdist(RFdist{roi1,roi2}*10,'Kernel');
         x = 0:2:100;
         ySix = pdf(pdSix,x);
         if roi1 ==roi2
            plot(x,ySix,'k-','LineWidth',2,'color',squeeze(Colors(roi1,3,:)));
         else
             plot(x,ySix,'k-','LineWidth',2,'color',[.4 .4 .4]);
         end
         
         YLIMs = ylim;
        text(60, YLIMs(2)*.9,['n=' num2str(numel(RFdist{roi1,roi2}))])
        if roi1==roi2
            title(ROI_names.(ROIs{roi2}))
        end
        if roi2==1
            ylabel(ROI_names.(ROIs{roi1}))
        end
    end
end

%% save figure

export_fig(FIG1,fullfile(Path,'RFCentroids'),'-pdf','-r200');
export_fig(FIG2,fullfile(Path,'RFAreas'),'-pdf','-r200');
export_fig(FIG3,fullfile(Path,'RFCentroidsDistance'),'-pdf','-r200');

end

function [Cent Pval Area Label] = extract_coords(unitid,probeid,sessionid,Probe_all)
    ind         = (Probe_all.Unit_ID==unitid) & (str2num(Probe_all.Probe_ID)==(probeid)) & (str2num(Probe_all.Session_ID)==(sessionid));
    [Pval cind] = min(Probe_all.RF_PValue(ind,:));
    Cent        = Probe_all.RF_Centroid(ind,:,cind);
    Area        = Probe_all.RF_Area(ind,cind)*100;
    Label       = Probe_all.Labels{ind};
end