clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
load('/Users/elhamb/switchdrive/EB/AllenBrainData/STOK_Average_iPDC_ff.98_MOrd10_Thalamus.mat');
%load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_Average_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs([1:7]));
%%
% Here, let's run a test :)

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
%PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);
PDC     = PDC(1:42,1:42,:,:);
for i = 1:size(PDC,1) % remove diagonal
    PDC(i,i,:,:) = 0;
end



%% overall output based on % output = intra/(intra+inter)
for roi= 1:7
    figure;
    fac = 20;
    ind1 = (roi-1)*6+1:roi*6;
    ind2 = setdiff(1:size(PDC,1),ind1);
    for i = 1:numel(ind1)

        subplot(6,3,(i-1)*3+2);
        imagesc(squeeze(nansum(PDC(ind2,ind1(i),:,Time>-.3),1))./(squeeze(nansum(PDC(:,ind1(i),:,Time>-.3),1))))
        axis xy
        caxis([-1.2 1.2])
        if i==1,title('inter');end

        subplot(6,3,(i-1)*3+1);
        imagesc(squeeze(nansum(PDC(ind1,ind1(i),:,Time>-.3),1))./(squeeze(nansum(PDC(:,ind1(i),:,Time>-.3),1))))
        axis xy
        caxis([-1.2 1.2])
        if i==1,title('intra');end
        ylabel(['L' num2str(i)],'fontweight','bold')

        subplot(6,3,(i-1)*3+3);
        Temp = squeeze(nansum(PDC(ind1,ind1(i),:,Time>-.3),1))./(squeeze(nansum(PDC(:,ind1(i),:,Time>-.3),1)));
        imagesc(Temp-mean(Temp(:,Time<0),2));
        %imagesc(+squeeze(nanmean(PDC(ind1,ind1(i),:,Time>-.3),1))-squeeze(nanmean(PDC(ind2,ind1(i),:,Time>-.3),1)));
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('intra % change');end
    end
    axes('position',[.45 .98 .1 .05]); axis off
    text(0,0,STOK_avg.ROIs{roi})
end
%% overall output
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

for roi= 1:7
    figure;
    fac = 2;
    ind1 = (roi-1)*6+1:roi*6;
    ind2 = setdiff(1:size(PDC,1),ind1);
    for i = 1:numel(ind1)

        subplot(6,3,(i-1)*3+2);
        imagesc(squeeze(nanmean(PDC(ind2,ind1(i),:,Time>-.3),1)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('inter');end

        subplot(6,3,(i-1)*3+1);
        imagesc(squeeze(nanmean(PDC(ind1,ind1(i),:,Time>-.3),1)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('intra');end
        ylabel(['L' num2str(i)],'fontweight','bold')

        subplot(6,3,(i-1)*3+3);
        imagesc(squeeze(nanmean(PDC(ind1,ind1(i),:,Time>-.3),1))-squeeze(nanmean(PDC(ind2,ind1(i),:,Time>-.3),1)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('intra-inter');end
    end
    axes('position',[.45 .98 .1 .05]); axis off
    text(0,0,STOK_avg.ROIs{roi})
end

%% overall input
for roi= 1:7
    figure;
    fac = 2;
    ind1 = (roi-1)*6+1:roi*6;
    ind2 = setdiff(1:size(PDC,1),ind1);
    for i = 1:numel(ind1)

        subplot(6,3,(i-1)*3+2);
        imagesc(squeeze(nanmean(PDC(ind1(i),ind2,:,Time>-.3),2)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('inter');end

        subplot(6,3,(i-1)*3+1);
        imagesc(squeeze(nanmean(PDC(ind1(i),ind1,:,Time>-.3),2)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('intra');end
        ylabel(['L' num2str(i)],'fontweight','bold')

        subplot(6,3,(i-1)*3+3);
        imagesc(squeeze(nanmean(PDC(ind1(i),ind1,:,Time>-.3),2))-squeeze(nanmean(PDC(ind1(i),ind2,:,Time>-.3),2)))
        axis xy
        caxis([-1.2 1.2]./fac)
        if i==1,title('intra-inter');end
    end
    axes('position',[.45 .98 .1 .05]); axis off
    text(0,0,STOK_avg.ROIs{roi})
end

