addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%%
load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','PARAFAC');

%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
addpath(genpath('../../ARC'));
IDs = fieldnames(StokALL);

% bootstrap
ROIsPerID = cellfun(@(x) StokALL.(x).ROIs,IDs,'uni',false);
nboots = 100;
bootsize = 5;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);

NComp = 3;

for b = 1:nboots
    clear Stok_sample;
    % make the sub-sample stok structure
    for s = 1:bootsize
        Stok_sample.(['S' num2str(s)]) = StokALL.(BootIDs{b,s});
    end

    disp(num2str(b))
    [PDC,~,Time,Freq,ROIs] = ExtractAllRoiPDC(Stok_sample,ROIs);
    %
    clear PDCFF PDCFB indFF indFB
    ind = 1;
    indl = reshape(1:36,[6 6]);
    for roi1 = 1:numel(ROIs)
        for roi2 = roi1+1:numel(ROIs)
            PDCFF{ind} = PDC{roi2,roi1};
            indFF(ind) = indl(roi2,roi1);
            PDCFB{ind} = PDC{roi1,roi2};
            indFB(ind) = indl(roi1,roi2);
            ind = ind+1;
        end
    end

    PDCFFtemp = cellfun(@(x) mean(x,5),PDCFF,'uni',false);
    PDCFBtemp = cellfun(@(x) mean(x,5),PDCFB,'uni',false);
    PDCPulled = cat(2,PDCFFtemp,PDCFBtemp);
    PDCPulled = cat(5,PDCPulled{:});
    indTotal = [indFF indFB];
    
    TW = Time>-0.5 & Time<1;
    [model_temp,it(b),err(b),corcondia(b)]=parafac(permute(PDCPulled(:,:,:,TW,:),[4 1 2 3 5]),NComp,[0 1 0 0 NaN],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    model{b} = model_temp([5 2 3 4 1]);
    temp_time = Time(TW);
end
save([FileName 'PARAFAC_covtemp_' num2str(NComp)],'model','NComp','temp_time','BootIDs','nboots','ROIs','IDs','Freq','indTotal','it','err','corcondia')
load([FileName 'PARAFAC_covtemp_' num2str(NComp)])

%% cluster the components  500 x 500 similarity
% for example here I calculate correlation for the first mode : this does
% not work

% We have different modes of data ...
for i = 1:nboots
    for j = 1:nboots
        [maps(i,j,:),sscore(i,j),Corrs(i,j,:),~,CorrPvals(i,j,:,:),Corrvals(i,j,:,:)] = PARAFACmodelComp(model{i},model{j},1:5);
    end
end

[~,RefInd] = max(mean(sscore)); % find the reference bootstrap
Corrsref = squeeze(Corrs(RefInd,:,:));

% reordering the models according to the reference component
for i = 1:nboots
    [~,sscore_reord(i),Corrs_reord(i,:),Model_reord{i}] = PARAFACmodelComp(model{RefInd},model{i});
end
%%
% cluster components
NClust = 6;

Corrvals_org = reshape(permute(Corrvals,[1 3 2 4]),size(Corrvals,1)*size(Corrvals,3),size(Corrvals,2)*size(Corrvals,4));
Dists = 1-Corrvals_org;
Z = linkage(Dists);
c = cophenet(Z,Dists);
I = inconsistent(Z);
T = cluster(Z,'maxclust',NClust);
TH = 10;

bnum = repmat(1:nboots,[NComp 1])'; bnum = bnum(:);
CRCN = repmat(corcondia,[NComp 1])'; CRCN = CRCN(:);
cnum = repmat(1:NComp,[nboots 1]);  cnum = cnum(:);

Mmodel = cat(1,model{:});
figure;
for clust = 1:NClust
    Boot_Clust = [bnum(T==clust & CRCN>TH) cnum(T==clust & CRCN>TH)];
    clear ttemp;
    for i = 1:size(Boot_Clust)
        temp = Mmodel(Boot_Clust(i,1),:);
        ttemp{i} = cellfun(@(x) x(:,Boot_Clust(i,2)),temp,'uni',false);
    end
    if exist('ttemp')
        ttemp = cat(1,ttemp{:});
        ttt= ttemp(:,5);
        %ttt = cellfun(@(x) (x-mean(x(temp_time<0 & temp_time>-.3)))./mean(x(temp_time<0 & temp_time>-.3)),ttt,'uni',false);
        %ttt = cellfun(@(x) x./norm(x),ttt,'uni',false);
        Temp = cat(2,ttt{:});%mean(cat(2,ttt{:}),2);
        subplot(NClust,1,clust),plot(Temp);
        title(['N = ' num2str(size(ttemp,1))])
    end
end

%legend({'C1','C2','C3','C4','C5'})

