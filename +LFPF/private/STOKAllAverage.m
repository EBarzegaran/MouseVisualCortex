function STOKAllAverage(STOKAll,ROIs_Select,savefig,Path,PDCMethod)

%% Organize ROIs and their order
Sessions_ID = fieldnames(STOKAll);

ROIs = cellfun(@(x) (STOKAll.(x).ROIs),Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(2,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

[C,ia,ib] = intersect(ROIs_Select,ROIs);
[~,I] = sort(ia);
ROIs = C(I);
ROI_counts = ROI_counts(ib(I));

%% Make the grand average matrix with all ROIS/Full model
Dims = size(STOKAll.(Sessions_ID{1}).PDC);
%
STOK_AV = zeros(6*numel(ROIs),6*numel(ROIs),Dims(3),Dims(4));
ROI_num = zeros(6*numel(ROIs),6*numel(ROIs));
for S = 1:numel(Sessions_ID)
   PDC = STOKAll.(Sessions_ID{S}).PDC;
   rois = STOKAll.(Sessions_ID{S}).ROIs;
   [C,ia] = intersect(ROIs_Select,rois);
   [~,I] = sort(ia);
   rois = C(I);
   Labels = ia(I);
   
   indices = arrayfun(@(x) [(x-1)*6+1:x*6],Labels,'uni',false);
   indices = cat(2,indices{:});
   
   STOK_AV(indices,indices,:,:) = STOK_AV(indices,indices,:,:) + PDC;
   ROI_num(indices,indices,:,:) = ROI_num(indices,indices,:,:) + 1;
   
end

STOK_AV = squeeze(STOK_AV)./ROI_num;
%--------------------------------------------------------------------------
saveresults = false;
if saveresults
    STOK_avg.PDC    = STOK_AV;
    STOK_avg.Time   = STOKAll.(Sessions_ID{1}).Times;
    STOK_avg.Freq   = STOKAll.(Sessions_ID{1}).f;
    STOK_avg.ROIs   = ROIs;
    STOK_avg.ROIs_count = ROI_counts;
    labels = cellfun(@(y) arrayfun(@(x) [y '_L' num2str(x)],1:6,'uni',false),ROIs,'uni',false);
    labels = cat(2,labels{:});
    STOK_avg.labels  = labels;
    save(fullfile(Path,'Fullmodel',['STOK_Average_' PDCMethod]),'STOK_avg');
    clear STOK_avg;
end
%% plot ROI/Layers connectivity matrix
Time = STOKAll.(Sessions_ID{1}).Times;
STOK_AVN = (STOK_AV - mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4))./mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4);

Data = squeeze(mean(STOK_AV(:,:,:,:),3));


A{1}      = squeeze(mean(Data(:,:,Time<.0 & Time>-.3),3));% prestimulus
A{2}      = squeeze(mean(Data(:,:,Time>.05 & Time<.25),3));% poststimulus
A{3}      = (A{2}-A{1})./A{1};% percent change %squeeze(mean(mean(STOK_AVN(:,:,:,Time>.05 & Time<.3),3),4));%

Titles = {'Prestimulus','Poststimulus','+ % Change','- %change'};
FIG = figure;
set(gcf,'unit','inch','color','w','position',[0 0 25 5.5])

for i = 1:3
    subplot(1,3,i),
    imagesc(A{i});
    %axis xy;
    if i==3
        M = max(abs(A{i}(:)));
        caxis([-M M]/1.5)
    else
        M = max(A{i}(:));
        m = min(A{i}(:));
        caxis([m M/1.5])
    end
    if i==1
        xlabel('Source')
        ylabel('Target')
    end
    hold on;

    set(gca,'xtick',3:6:numel(ROIs)*6,'xticklabel',ROIs,'ytick',3:6:numel(ROIs)*6,'yticklabel',ROIs)
    for r = 0:numel(ROIs)-1
        vline(r*6+.5,'w-')
        hline(r*6+.5,'w-')
    end
    colorbar
    title(Titles{i})
end

if savefig
    export_fig(FIG,fullfile(Path,['Average_All_STOK_Transient_' PDCMethod]),'-pdf');
end

%% plot the graphs: within ROI, laminar connectivity

LLabels = arrayfun(@(x) [ num2str(x)],1:6,'uni',false);
Colors  = brewermap(6,'Spectral');

AN = cellfun(@(x) x./max(x(:)),A,'uni',false);
%AN{4} = -1*AN{3};

FIG = figure;
set(gcf,'unit','inch','color','w','position',[1 1 numel(AN)*2+2 25])
for roi = 1:numel(ROIs)
    ind = (roi-1)*6+1:roi*6;
    for i = 1:numel(AN)
        %-----------------------------
        subplot(numel(ROIs),numel(AN),(roi-1)*numel(AN)+i);
        a = AN{i}(ind,ind);
        a(a<0)=0;
        a = a./max(a(:));
        plot_graph(a,LLabels,Colors,'down',.6)
        %------------------------------
        if roi ==1
            title(Titles{i});
        end
        if i==1
            ylabel(ROIs{roi});
        end
        
        if roi ==numel(ROIs)
            xlabel('Upward     Downward')
        end
    end
   
end

if savefig
    export_fig(FIG,fullfile(Path,['Average_All_Laminar_Transient_' PDCMethod]),'-pdf');
end

%% plot network structure
FIG = figure;
set(gcf,'unit','inch','color','w','position',[1 1 15 5])
Colors  = brewermap(8,'Spectral');

for i = 1:3
    subplot(1,3,i),
    for r1 = 1:numel(ROIs)
        for r2 = 1:numel(ROIs)
            a(r1,r2) = max(max(AN{i}((r1-1)*6+1:r1*6,(r2-1)*6+1:r2*6)));
        end
    end
%     a = squeeze(mean(mean(reshape(A{i},6,numel(ROIs),6,numel(ROIs)),1),3));
    a(a<0)=0;
    a = a./max(a(:));
    plot_graph(a,ROIs,Colors,'down',.5)
    title(Titles{i});
    xlabel('Feedback         Feedforward')
end

if savefig
    export_fig(FIG,fullfile(Path,['Average_All_ROIs_Transient_' PDCMethod]),'-pdf');
end
%%
Freq = STOKAll.(Sessions_ID{1}).f;

for roi = 1:numel(ROIs)
    ind = (roi-1)*6+1:roi*6;
   

    FIG = dynet_connplot(STOK_AVN(ind,ind,:,Time>-.1 & Time<.2),Time(Time>-.1 & Time<.2),Freq);
    set(FIG,'color','w','unit','inch','position',[0 0 25 15]);
    if savefig
        export_fig(FIG,fullfile(Path,'Fullmodel',['STOK_individual_Transient_' ROIs{roi} '_' PDCMethod]),'-pdf');
    end
    close all
end

end


