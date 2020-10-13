function PARAFAC_FC(StokALL,NComp,BootIDs,nboots,ROIs,DataPath,FigPath,FileName,redo)


if ~exist([FileName 'PARAFAC_covtemp_' num2str(NComp) '.mat'],'file') || redo
    for b = 1:nboots
            clear Stok_sample;
            % make the sub-sample stok structure
            for s = 1:size(BootIDs,2)
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
            %-------------------unfolded matrix only in between-areal connectivity dimesions-------------------------
             [model_temp,it(b),err(b),corcondia(b)]=parafac(permute((PDCPulled(:,:,:,TW,:)),[4 1 2 3 5]),NComp,[1e-7 10 0 0 NaN],repmat(2,1,5));% dimensions: connections x in x out x freq x time
            
             model{b} = model_temp([5 2 3 4 1]);
            %-------------------unfold the matrix also in laminar connectivity dimesions-------------------------
            %%%%%%% NOTE: It does not help to improve corcondia scores%%%%%
            % PDCP = reshape(PDCPulled(:,:,:,:,:),6*6,size(PDCPulled,3),size(PDCPulled,4),size(PDCPulled,5));
            % [model_temp,it(b),err(b),corcondia(b)]=parafac(permute((PDCP(:,:,TW,:)),[3 1 2 4]),NComp,[1e-7 10 0 0 NaN],repmat(2,1,5));
            
            % model{b} = model_temp([4 2 3 1]);  
            %-------------------------average or selecting laminar layers-------------------------
            
            %-------------------------TUCKER---------------------------
            % Factors,G,ExplX,Xm]=tucker(permute((PDCP(:,:,TW,:)),[3 1 2 4]),[4 8 3 4],[1e-7 10 0 0 NaN]);


            
            temp_time = Time(TW);
    end
    save(fullfile(DataPath,[FileName 'PARAFAC_covtemp_' num2str(NComp)]),'model','NComp','temp_time','BootIDs','nboots','ROIs','IDs','Freq','indTotal','it','err','corcondia')
else
    load(fullfile(DataPath,[FileName 'PARAFAC_covtemp_' num2str(NComp)]));
end
%% cluster the components  500 x 500 similarity
% for example here I calculate correlation for the first mode : this does
% not work
% We have different modes of data ...
if ~exist([FileName 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar.mat'],'file')
    for i = 1:nboots
        i
        for j = 1:nboots
            [maps(i,j,:),sscore(i,j),Corrs(i,j,:),~,CorrPvals(i,j,:,:),Corrvals(i,j,:,:)] = PARAFACmodelComp(model{i},model{j});
        end
    end

    [~,RefInd] = max(mean(sscore)); % find the reference bootstrap
    Corrsref = squeeze(Corrs(RefInd,:,:));

    % reordering the models according to the reference component
    for i = 1:nboots
        [~,sscore_reord(i),Corrs_reord(i,:),Model_reord{i}] = PARAFACmodelComp(model{RefInd},model{i});
    end
else
    load([FileName 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar.mat']);
end
%% how much of variance are explained by each component
%[1 4 3 2 5];

clear Model_vars;
for y = 1:NComp
    Model_vars{y} = cellfun(@(x) norm(x{5}(:,y)),Model_reord);
end



Model_vars=cat(1,Model_vars{:})';
Model_var_percent = Model_vars./sum(Model_vars,2);

[~,Comp_ord] = sort(mean(Model_var_percent),'descend');%[1 3 5 2 4];

FIG = figure;
set(FIG,'unit','inch','position',[0 0 2*NComp 4],'color','w')
bar(1:NComp,mean(Model_var_percent(:,Comp_ord)),.5);hold on;
errorbar(1:NComp,mean(Model_var_percent(:,Comp_ord)),var(Model_var_percent(:,Comp_ord)),'.','color','k')
xlabel('Component');
ylabel('Variance explained')
xlim([.5 NComp+.5])
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_VarianceExplained']),'-pdf','-r200')

%% if you want to check the components' Consistency, you should use the non-ordered ones: corrs and maps

% (1) Everything according to the reference bootstrap
for i = 1:nboots
    for j = 1:nboots
        map1 = squeeze(maps(RefInd,i,:));
        map2 = squeeze(maps(i,j,:));
        Corrs_reord2(i,j,:) = Corrs(i,j,map1);% or squeeze(maps(RefInd,j,:))
        CorrPvals_reord2(i,j,:,:) = CorrPvals(i,j,map1,:);% or squeeze(maps(RefInd,j,:))
    end
end

Comp_PerSig = squeeze(sum(reshape(CorrPvals_reord2,nboots*nboots,NComp,5)<.05)/(nboots.^2));


FIG = figure;
set(FIG,'unit','inch','position',[0 0 2*NComp 4],'color','w')
bar(1:NComp,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),.5);hold on;
errorbar(1:NComp,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),var(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),'.','color','k')
xlim([.5 NComp+.5])
ylim([0.0 1.2])
xlabel('Component');
ylabel('Average inter boostrap Correlations')
for c = 1:NComp
    text(c,1,['Consistency=' num2str(round(mean(Comp_PerSig(Comp_ord(c),:)),2))],'HorizontalAlignment','center')
end
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'Consistency']),'-pdf','-r200')
save(fullfile(DataPath,[FileName 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar']),'maps','Corrs','CorrPvals','sscore','Model_reord','Comp_ord','Comp_PerSig','RefInd')

end