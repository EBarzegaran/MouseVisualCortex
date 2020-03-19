function STOKAllAverage(STOKAll,ROIs_Select,savefig,Path,SaveName)
% STOKALL is a structure with animal ID as field names and each
% STOKALL.Session_ID is a structure itself with the following fields: PDC, f,
% Times, ROIs

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
ROIsizes = contains(ROIs_Select,'VIS')*3+3;
ROIindices = [0 cumsum(ROIsizes)];
%
STOK_AV = zeros(sum(ROIsizes),sum(ROIsizes),Dims(3),Dims(4));
ROI_num = zeros(sum(ROIsizes),sum(ROIsizes));
for S = 1:numel(Sessions_ID)
   PDC = STOKAll.(Sessions_ID{S}).PDC;
   rois = STOKAll.(Sessions_ID{S}).ROIs;
   [C,ia] = intersect(ROIs_Select,rois);
   [~,I] = sort(ia);
   rois = C(I);
   Labels = ia(I);
   
   %indices = arrayfun(@(x) [(x-1)*6+1:x*6],Labels,'uni',false);
   %indices = cat(2,indices{:});
   indices  =  arrayfun(@(x) ROIindices(x)+1:ROIindices(x+1),Labels,'uni',false);
   indices = cat(2,indices{:});
   
   STOK_AV(indices,indices,:,:) = STOK_AV(indices,indices,:,:) + PDC(:,:,:,1:Dims(4));
   ROI_num(indices,indices,:,:) = ROI_num(indices,indices,:,:) + 1;
   
end

STOK_AV = squeeze(STOK_AV)./ROI_num;
%--------------------------------------------------------------------------
saveresults = true;
if saveresults
    STOK_avg.PDC    = STOK_AV;
    STOK_avg.Time   = STOKAll.(Sessions_ID{1}).Times;
    STOK_avg.Freq   = STOKAll.(Sessions_ID{1}).f;
    STOK_avg.ROIs   = ROIs_Select;
    STOK_avg.ROIs_count = ROI_counts;
    labels = arrayfun(@(y) arrayfun(@(x) [ROIs_Select{y} '_L' num2str(x)],1:ROIsizes(y),'uni',false),1:numel(ROIs_Select),'uni',false);
    labels = cat(2,labels{:});
    STOK_avg.labels  = labels;
    save(fullfile(Path,'Fullmodel',['STOK_Average_' SaveName]),'STOK_avg');
    clear STOK_avg;
end
%% plot ROI/Layers connectivity matrix
if false
    Time = STOKAll.(Sessions_ID{1}).Times;
    STOK_AVN = (STOK_AV - mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4))./mean(STOK_AV(:,:,:,Time<0 & Time>-.3),4);

    Data = squeeze(mean(STOK_AV(:,:,:,:),3));


    A{1}      = squeeze(mean(Data(:,:,Time<.0 & Time>-.3),3));% prestimulus
    A{2}      = squeeze(mean(Data(:,:,Time>.05 & Time<.25),3));% poststimulus
    A{3}      = (A{2}-A{1})./A{1};% percent change %squeeze(mean(mean(STOK_AVN(:,:,:,Time>.05 & Time<.3),3),4));%

    Titles = {'Prestimulus','Poststimulus','+ % Change','- %change'};
    FIG = figure;
    set(gcf,'unit','inch','color','w','position',[0 0 25 5.5])
    ROItick = (ROIindices(1:end-1)+ROIindices(2:end))/2;
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

        set(gca,'xtick',ROItick,'xticklabel',ROIs_Select,'ytick',ROItick,'yticklabel',ROIs_Select)
        xtickangle(45)
        for r = 2:numel(ROIs_Select)
            vline(ROIindices(r)+.5,'w-')
            hline(ROIindices(r)+.5,'w-')
        end
        colorbar
        title(Titles{i})
    end

    if savefig
        export_fig(FIG,fullfile(Path,['Average_All_STOK' SaveName]),'-pdf');
    end

    %% plot the graphs: within ROI, laminar connectivity

    %LLabels = arrayfun(@(x) [ num2str(x)],1:6,'uni',false);
    Colors  = brewermap(6,'Spectral');

    AN = cellfun(@(x) x./max(x(:)),A,'uni',false);
    %AN{4} = -1*AN{3};

    FIG = figure;
    set(gcf,'unit','inch','color','w','position',[1 1 numel(AN)*2+2 25])
    for roi = 1:numel(ROIs_Select)
        ind = ROIindices(roi)+1:ROIindices(roi+1);
        for i = 1:numel(AN)
            %-----------------------------
            subplot(numel(ROIs_Select),numel(AN),(roi-1)*numel(AN)+i);
            a = AN{i}(ind,ind);
            a(a<0)=0;
            a = a./max(a(:));
            LLabels = arrayfun(@(x) [ num2str(x)],1:ROIsizes(roi),'uni',false);
            plot_graph(a,LLabels,Colors,'down',.6)
            %------------------------------
            if roi ==1
                title(Titles{i});
            end
            if i==1
                ylabel(ROIs_Select{roi});
            end

            if roi ==numel(ROIs_Select)
                xlabel('Upward     Downward')
            end
        end

    end

    if savefig
        export_fig(FIG,fullfile(Path,['Average_All_Laminar' SaveName]),'-pdf');
    end

    %% plot network structure
    FIG = figure;
    set(gcf,'unit','inch','color','w','position',[1 1 15 5])
    Colors  = brewermap(11,'Spectral');

    for i = 1:3
        subplot(1,3,i),
        for r1 = 1:numel(ROIs_Select)
            for r2 = 1:numel(ROIs_Select)
                a(r1,r2) = max(max(AN{i}(ROIindices(r1)+1:ROIindices(r1+1),ROIindices(r2)+1:ROIindices(r2+1))));
            end
        end
    %     a = squeeze(mean(mean(reshape(A{i},6,numel(ROIs),6,numel(ROIs)),1),3));
        a(a<0)=0;
        a = a./max(a(:));
        plot_graph(a,ROIs_Select,Colors,'down',.7)
        title(Titles{i});
        xlabel('Feedback         Feedforward')
    end

    if savefig
        export_fig(FIG,fullfile(Path,['Average_All_ROIs' SaveName]),'-pdf');
    end
    %%
    Freq = STOKAll.(Sessions_ID{1}).f;

    for roi = 1:numel(ROIs_Select)
        ind = ROIindices(roi)+1:ROIindices(roi+1);


        FIG = dynet_connplot(STOK_AVN(ind,ind,:,Time>-.1 & Time<.2),Time(Time>-.1 & Time<.2),Freq);
        set(FIG,'color','w','unit','inch','position',[0 0 25 15]);
        if savefig
            export_fig(FIG,fullfile(Path,'Fullmodel',['STOK_individual_' ROIs_Select{roi} SaveName]),'-pdf');
        end
        close all
    end
end
end


