function STOKAllLayerStats(STOKAll,ROIs_Select,savefig,Path,PDCMethod)

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

Time = STOKAll.(Sessions_ID{1}).Times;
STOK_AVN = (STOK_AV - mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4))./mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4);

Data = squeeze(mean(STOK_AV(:,:,:,:),3));


A{1}      = squeeze(mean(Data(:,:,Time<.0 & Time>-.3),3));% prestimulus
A{2}      = squeeze(mean(Data(:,:,Time>.05 & Time<.25),3));% poststimulus
A{3}      = (A{2}-A{1})./A{1};% percent change %squeeze(mean(mean(STOK_AVN(:,:,:,Time>.05 & Time<.3),3),4));%

%%
AT = reshape(A{3},6,numel(ROIs),6,numel(ROIs));
ind = 1;
for i = 1:numel(ROIs)
    for j = 1:numel(ROIs)
        if i ~= j 
            for k = 1:6
                if ~isnan(sum(AT(:,i,k,j)))
                    LamD(:,ind) = AT(:,i,k,j);
                    InfoStrct.Layer(ind) = k;
                    InfoStrct.Source(ind) = j;
                    InfoStrct.Target(ind) = i;
                    ATR(:,k,ind)=AT(:,i,k,j);
                    ind = ind+1;
                    
                end
            end
        end
    end
end

%% First Figure
FIG = figure;
set(gcf,'unit','inch','color','w','position',[5 5 8 7]);

imagesc(mean(ATR,3));
colormap(jmaColors('coolhot'));
caxis([0 max(max(mean(ATR,3)))]);
set(gca,'xtick',1:6,'ytick',1:6);
xlabel('Source');
ylabel('Target');
title('Between Area Laminar Connectivity')
colorbar;

if savefig
    export_fig(FIG,fullfile(Path,'Fullmodel',['Layers_LaminarConnectivityAverage_' PDCMethod]),'-pdf');
end
%% Second figure

AT = LamD;
AT = (AT-min(AT))./(max(AT)-min(AT));


FIG = figure;
set(gcf,'unit','inch','color','w','position',[0 1 8 10]);

eucD = pdist(AT','euclidean');
clustTreeEuc = linkage(eucD,'average');
cophenet(clustTreeEuc,eucD)
% figure,
% [h,nodes] = dendrogram(clustTreeEuc,7);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = [];
hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',1.0); %For normalized
%hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',.05);
%---------------------Sort the clusters---------------------------
ATM = [];
MI = [];
for i = 1:max(hidx)
    ATM(:,i) = median(LamD(:,hidx==i),2);
    [~,MI(i)]= max(ATM(:,i)); 
end
[~,SI] = sort(MI);
%SI = 1:numel(MI);
ATM = ATM(:,SI);
%---------------------Organize the input-----------------------------
ind = 1;
IL = [];
for i = 1:max(hidx)
    Ind = find(hidx==SI(i));
    AT2(:,ind:ind+numel(Ind)-1) =  AT(:,Ind);
    ind = ind+numel(Ind);
    IL(i)=ind;
end
subplot(2,1,1);imagesc(AT2)
%colormap('hot')
ILS = [0 IL];
for i = 1:numel(IL-1)
    line([IL(i) IL(i)]-.5,[0 6.5],'color','y','linewidth',2);
    XTL(i) = round((ILS(i)+ILS(i+1))/2);
end
set(gca,'xtick',XTL,'xticklabel',arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false));
colorbar;
%--------------------------------------------------------------------
subplot(2,1,2)
Colors  = brewermap(max(hidx),'Spectral');
% for i = 1:max(hidx)
%     [~,m] = max(ATM(:,i));
%     plot(ATM(:,i)+(.3*i),'--o','linewidth',1.5,'color',Colors(i,:));hold on
% end
% legend(arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false))
% set(gca,'xtick',1:6,'xticklabel',arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
% xlim([.5 7])
imagesc(ATM);
colormap(jmaColors('coolhot')); 
colorbar;
for i = 1:numel(IL-1)
    line([i i]-.5,[0 6.5],'color','y','linewidth',2);
end
set(gca,'xtick',1:max(hidx),'xticklabel',arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false));

if savefig
    export_fig(FIG,fullfile(Path,'Fullmodel',['Layers_Clustering_' PDCMethod]),'-pdf');
end

%% --------------------------------------------------------------------------
% Third figure
clear Source Target Layer
for h = 1:max(hidx)
    [uv,~,idx] = unique(InfoStrct.Source(hidx==h));
    n = accumarray(idx(:),1);
    Source(h,uv)=n;
    
    [uv,~,idx] = unique(InfoStrct.Target(hidx==h));
    n = accumarray(idx(:),1);
    Target(h,uv)=n;
    
    [uv,~,idx] = unique(InfoStrct.Layer(hidx==h));
    n = accumarray(idx(:),1);
    Layer(h,uv)=n;
    
end
%-----------------------------------------------------------------
Colors = brewermap(max(hidx),'Spectral');

FIG = figure;
set(gcf,'unit','inch','color','w','position',[8 1 8 15]);
Clabel = MI(SI);
subplot(3,1,1)
B = bar(Layer(SI,:)','stacked','FaceColor','flat');
for k = 1:size(Layer,1)
    B(k).CData = Colors((k),:);
end
legend(arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false))
set(gca,'xticklabel',arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
xlabel('Source Layer')
ylabel('#')


subplot(3,1,2)
B = bar(Source(SI,:)','stacked','FaceColor','flat');
for k = 1:size(Layer,1)
    B(k).CData = Colors((k),:);
end
set(gca,'xticklabel',ROIs)
xlabel('Source ROI')
ylabel('#')

subplot(3,1,3)
B = bar(Target(SI,:)','stacked','FaceColor','flat');
for k = 1:size(Layer,1)
    B(k).CData = Colors((k),:);
end
set(gca,'xticklabel',ROIs)
xlabel('Target ROI')
ylabel('#')

if savefig
    export_fig(FIG,fullfile(Path,'Fullmodel',['Layers_ClusterHist1_' PDCMethod]),'-pdf');
end

%% -----------------------------------------------------------------------
% Fourth figure

FIG = figure;
set(gcf,'unit','inch','color','w','position',[16 1 8 15]);
Colors = brewermap(6,'Spectral');

subplot(3,1,1)
B = bar(Layer(SI,:),'stacked','FaceColor','flat');
for k = 1:size(Layer,2)
    B(k).CData = Colors(k,:);
end
legend(arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
set(gca,'xticklabel',arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false))
title('Source Layer')
ylabel('#')

Colors = brewermap(8,'Spectral');
subplot(3,1,2)
B = bar(Source(SI,:),'stacked','FaceColor','flat');
for k = 1:size(Source,2)
    B(k).CData = Colors(k,:);
end
set(gca,'xticklabel',arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false))
legend(ROIs)
title('Source ROI')
ylabel('#')

subplot(3,1,3)
B = bar(Target(SI,:),'stacked','FaceColor','flat');
for k = 1:size(Target,2)
    B(k).CData = Colors(k,:);
end
set(gca,'xticklabel',arrayfun(@(x) ['C' num2str(x)],1:max(hidx),'uni',false))
legend(ROIs)
title('Target ROI')

ylabel('#')
if savefig
    export_fig(FIG,fullfile(Path,'Fullmodel',['Layers_ClusterHist2_' PDCMethod]),'-pdf');
end

%% Fifth figure
FIG = figure;
set(gcf,'unit','inch','color','w','position',[5 5 8 7]);


for i = 1:numel(ROIs)
    for j = 1:numel(ROIs)
        if i ~= j 
            for k = 1:6
                    
                ind = find((InfoStrct.Layer == k) & (InfoStrct.Source==j) & (InfoStrct.Target==i));
                if ~isempty(ind)
                    ClustM(i,(j-1)*6+k) = find(SI == hidx(ind));
                end
            end
        end
    end
end


imagesc(ClustM)
for j = 1:numel(ROIs)
        line([(j-1)*6 (j-1)*6]+.5,[0 9],'color','k','linewidth',2);
        line([0 numel(ROIs)*6],[j j]+.5,'color','k','linewidth',2);
end
Colors = brewermap(max(hidx),'Spectral');
Colors = [0.5 0.5 0.5; Colors];
colormap(Colors);
set(gca,'ytick',(1:numel(ROIs)),'yticklabels',ROIs,'xtick',(3.5:6:numel(ROIs)*6),'xticklabels',ROIs);
cbh = colorbar;
set(cbh,'XTick',(0:(max(hidx)-1)/max(hidx):max(hidx))+.45,'XTickLabel',arrayfun(@(x) ['C' num2str(x)],0:max(hidx),'uni',false));
title('Cluster distributions')
xlabel('Source');
ylabel('Target');
if savefig
    export_fig(FIG,fullfile(Path,'Fullmodel',['Layers_Cluster distributions_' PDCMethod]),'-pdf');
end
end


