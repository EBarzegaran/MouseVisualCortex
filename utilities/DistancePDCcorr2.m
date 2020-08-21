function Results = DistancePDCcorr2(StokALL,P_th,Path,CondName)

    %P_th        = .2; 
    IDs         = fieldnames(StokALL);
    PDC_All     = [];
    Dist_All    = [];
    Pval_All    = [];
    Laminar_All    = [];
    Diag_All    = [];
    Col_All     = [];
    ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
    load ROInames;
    ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
    NROIs = numel(ROIs);
    Dists = zeros(NROIs,NROIs,6,numel(IDs)-1);
    RFDists = zeros(NROIs,NROIs,6,numel(IDs)-1);
    PDC_All = zeros(NROIs,NROIs,6,numel(IDs)-1,100,751);
    for id = 1:numel(IDs)-1
        PDC_temp    = StokALL.(IDs{id}).PDC;
        
        nrois = numel(StokALL.(IDs{id}).ROIs);
        roiind = find(ismember(ROIs,StokALL.(IDs{id}).ROIs));
        for roi1 = 1:nrois
            IND1 = (roi1-1)*6+1:roi1*6;
            for roi2 = roi1+1:nrois
                IND2 = (roi2-1)*6+1:roi2*6;
                Dists(roiind(roi1),roiind(roi2),:,id) = diag(StokALL.(IDs{id}).ProbeDist(IND1,IND2));
                RFDists(roiind(roi1),roiind(roi2),:,id) = diag(StokALL.(IDs{id}).ProbeRFDist(IND1,IND2));
                pdc_temp = reshape(PDC_temp(IND1,IND2,:,:),36,size(PDC_temp,3),size(PDC_temp,4));
                PDC_All(roiind(roi1),roiind(roi2),:,id,:,:) = pdc_temp(1:6:end,:,:);
            end
        end
        
    end
    Dists = Dists+permute(Dists,[2 1 3 4]);
    RFDists = RFDists+permute(RFDists,[2 1 3 4]);
    PDC_All = PDC_All + permute(PDC_All,[2 1 3 4 5 6]);
    %% physical distance
    figure,
    Dists(Dists==0)=NaN;
    RFDists(RFDists==0)=NaN;
    PDC_All(PDC_All==0) = NaN;
    subplot(1,2,1),imagesc(nanmean(nanmean(Dists(:,:,:,:),3),4))
    caxis([200 3000])
    subplot(1,2,2),imagesc(nanmean(nanmean(RFDists(:,:,:,:),3),4))
    %colormap()
    
    
    set(gca,'xtick',1:6,'xticklabels',ROISN,'ytick',1:6,'yticklabels',ROISN)
    %% plot V1 effect change
    PDC_All = (PDC_All- mean(PDC_All(:,:,:,:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.5),6))./mean(PDC_All(:,:,:,:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.5),6);
    %
    for f = 1:100
        for t = 100:751
            pdc = PDC_All(:,:,:,:,f,t);
            Dists = Dists(:);
            RFDists = RFDists(:);
            pdc = pdc(:);
           [C_Dists(f,t),P_Dists(f,t)] = corr(Dists(~isnan(Dists)),pdc(~isnan(pdc)));
           [C_RFDists(f,t),P_RFDists(f,t)] = corr(RFDists(~isnan(RFDists)&~isnan(pdc)),pdc(~isnan(pdc)&~isnan(RFDists)));
        end
    end
    
    %%
    FIG = figure;
    subplot(2,1,1),
    IM = imagesc(StokALL.S766640955.Times,StokALL.S766640955.f,C_Dists);
    set(IM,'alphadata',(P_Dists<0.025)*.8+.2);
    axis xy;
    xlim([-.5 1])
    caxis([-.2 .2])
    colormap(jet)
    line([0 0],[1 100],'color','r','linestyle','--','linewidth',1.5)
    title('Correlation with Distance');ylabel('Frequency(Hz)');
    colorbar
    
    subplot(2,1,2),
    IM = imagesc(StokALL.S766640955.Times,StokALL.S766640955.f,C_RFDists);
    set(IM,'alphadata',(P_RFDists<0.025)*.8+.2);
    axis xy;
    xlim([-.5 1])
    caxis([-.2 .2])
    colormap(jet)
    line([0 0],[1 100],'color','r','linestyle','--','linewidth',1.5)
    title('Correlation with RF Distance');ylabel('Frequency(Hz)'); xlabel('Time(S)')
    colorbar;
    
    %export_fig(FIG,fullfile(Path,['Distance_PDC_Corr_original']),'-pdf','-r200');
    export_fig(FIG,fullfile(Path,['Distance_PDC_Corr_evoked']),'-pdf','-r200');
end