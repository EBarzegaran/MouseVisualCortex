function STOKAllDirectStats(STOKAll,ROIs_Select,savefig,Path,PDCMethod)

%% Organize ROIs and their order
Sessions_ID = fieldnames(STOKAll);

ROIs = cellfun(@(x) (STOKAll.(x).ROIs),Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(2,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

[C,ia,ib] = intersect(ROIs_Select,ROIs);
[~,I] = sort(ia);
ROIs = C(I);
ROI_counts = ROI_counts(ib(I));

%%
Dims = size(STOKAll.(Sessions_ID{1}).PDC);
% Make the grand matrix with all ROIS
STOK_A = zeros(numel(Sessions_ID),6*numel(ROIs),6*numel(ROIs),Dims(3),Dims(4));
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
   
   STOK_A(S,indices,indices,:,:) = PDC;%STOK_AV(indices,indices,:,:) + PDC;
   ROI_num(indices,indices,:,:) = ROI_num(indices,indices,:,:) + 1;
   
end

%%
Time = STOKAll.(Sessions_ID{1}).Times;
%STOK_AN = (STOK_A - mean(STOK_A(:,:,:,:,Time<0 & Time>-.3),5))./mean(STOK_A(:,:,:,:,Time<0 & Time>-.3),5);

STOK_AM = squeeze(mean(mean(reshape(STOK_A,size(STOK_A,1),6,numel(ROIs),6,numel(ROIs),size(STOK_A,4),size(STOK_A,5)),2),4));
ROI_numM = squeeze(mean(mean(reshape(ROI_num,6,numel(ROIs),6,numel(ROIs)),1),3));
NBoots = 200;
% 1000 bootstrap over animals
for r1 = 1:numel(ROIs)
    for r2= 1:numel(ROIs)
        r2
        %tic
        if ROI_numM(r1,r2)>1
            Boots = randi(ROI_numM(r1,r2),NBoots,ROI_numM(r1,r2));
            for B = 1:NBoots
                STOK_AMB(B,r1,r2,:,:) = mean(STOK_AM(Boots(B,:),r1,r2,:,:),1);
            end
        end
        %toc
    end
end

STOK_B = squeeze(mean(STOK_AMB));

TWIN = find((Time>.1 & Time<.2));%numel(Time);%200:350;
dynet_connplot(STOK_B(:,:,:,TWIN),Time(TWIN),5:101,ROIs,[],[],[],1);
STOK_BN = (STOK_B - mean(STOK_B(:,:,:,Time<0 & Time>-.5),4))./mean(STOK_B(:,:,:,Time<0 & Time>-.3),4);
dynet_connplot(STOK_BN(:,:,:,TWIN),Time(TWIN),5:101,ROIs,[],[],[],1);

STOK_BN(STOK_BN<0)=0;
STOK_B_d = STOK_BN - permute(STOK_BN,[2 1 3 4]);
dynet_connplot(STOK_B_d(:,:,:,TWIN),Time(TWIN),5:101,ROIs);
colormap(jmaColors('coolhotcortex'))
% dynet_connplot(STOK_AVNM(:,:,:,Time>-.1 & Time<.3),Time(Time>-.1 & Time<.3),5:100,ROIs);
% dynet_connplot(STOK_AVNM_d(:,:,:,Time>-.1 & Time<.3),Time(Time>-.1 & Time<.3),5:100,ROIs);

a = squeeze(nanmean(nanmean(STOK_BN(:,:,:,Time>.1 & Time<.2),4),3));
Colors  = brewermap(8,'Spectral');
figure,plot_graph(a,ROIs,Colors,'down',.5)
%%
 STOK_BNM = mean(STOK_BN(:,:,:,TWIN),4);
 figure,
for i = 1:8
for j = i:8
subplot(8,8,(i-1)*8+j),plot(squeeze(STOK_BNM(i,j,:)),'b','linewidth',1.5)
hold on;plot(squeeze(STOK_BNM(j,i,:)),'r','linewidth',1.5)
ylim([-.3 1])
if i==1, title(ROIs{j}); end
end
end

%% distribution plots: directionality, up/downward - feedback/feedforward

% Data = squeeze(mean(STOK_A,4));
% STOK_AN = (STOK_A - mean(STOK_A(:,:,:,:,Time<0 & Time>-.3),5))./mean(STOK_A(:,:,:,:,Time<0 & Time>-.3),5);
% 
% A{1}      = squeeze(mean(Data(:,:,:,Time<.0 & Time>-.3),4));% prestimulus
% A{2}      = squeeze(mean(Data(:,:,:,Time>.05 & Time<.3),4));% poststimulus
% A{3}      = squeeze(mean(mean(STOK_AN(:,:,:,:,Time>.05 & Time<.3),4),5));%(A{2}-A{1})./A{1};% percent change
% 
% FIG  = figure;
% set(gcf,'unit','inch','color','w','position',[1 1 15 5])
% 
% for i = 1:3
%     subplot(1,3,i)
%     for roi = 1:numel(ROIs)
%         ind = (roi-1)*6+1:roi*6;
%         a = A{i}(:,ind,ind);
%         a = reshape(a,size(a,1),6*6);
%         I1 = tril(ones(6),-1);% upward
%         I2 = triu(ones(6),1);% downward
%         b1 = arrayfun(@(x) a(x,I1(:)==1),1:size(a,1),'uni',false);b1 = cat(2,b1{:});b1(b1==0)=[];
% 
%         b2 = arrayfun(@(x) a(x,I2(:)==1),1:size(a,1),'uni',false);b2 = cat(2,b2{:});b2(b2==0)=[];
%         M(roi,1)= nanmean(b1);
%         M(roi,2)= nanmean(b2);
%         SEM(roi,1)= nanstd(b1)./sqrt(numel(b1));
%         SEM(roi,2)= nanstd(b2)./sqrt(numel(b2));
%     end
% 
%     B = bar(M);hold on;
%     %errorbar(M(:,1),SEM(:,1))
%     if i==3, legend('Upward','Downward');end
%     set(gca,'xticklabels',ROIs);
%     xtickangle(45);
%     if i <3
%         ylim([0 .08])
%     end
%     title(Titles{i});
% end
% 
% if savefig
%     export_fig(FIG,fullfile(Path,['Average_All_ROIs_Transient_' PDCMethod '_directionality']),'-pdf');
% end
dynet_connplot
end