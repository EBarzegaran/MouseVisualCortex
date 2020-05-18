function Results = DistancePDCcorr(StokALL,P_th,Path,CondName)

    %P_th        = .2; 
    IDs         = fieldnames(StokALL);
    PDC_All     = [];
    Dist_All    = [];
    Pval_All    = [];
    Laminar_All    = [];
    Diag_All    = [];
    Col_All     = [];
    
    for id = 1:numel(IDs)
        id
        inter_connect = kron(eye(numel(StokALL.(IDs{id}).ROIs)),ones(6));
        %
        PDC_temp    = StokALL.(IDs{id}).PDC;
        NNode       = size(PDC_temp,1);
        PDC_temp    = reshape(PDC_temp,[NNode*NNode size(PDC_temp,3) size(PDC_temp,4)]);
        Dist_temp   = reshape(StokALL.(IDs{id}).ProbeDist,NNode*NNode,1);
        Pval_temp   = reshape((StokALL.(IDs{id}).ProbeRFPval<P_th).*(StokALL.(IDs{id}).ProbeRFPval'<P_th),NNode*NNode,1);
        Laminar_temp = reshape(inter_connect,NNode*NNode,1);
        Diag_temp   = reshape(eye(NNode),NNode*NNode,1);
        %Col_temp    = zeros(NNode);
        %Col_temp(7:end,1:6) = 1;
        %Col_temp    = tril(1-inter_connect); % FF
        Col_temp    = triu(1-inter_connect); % FB
        Col_temp    = reshape(Col_temp,NNode*NNode,1); 
        
        PDC_All     = cat(1,PDC_All,PDC_temp);
        Dist_All    = cat(1,Dist_All,Dist_temp);
        Diag_All    = cat(1,Diag_All,Diag_temp);
        Laminar_All = cat(1,Laminar_All,Laminar_temp);
        Col_All     = cat(1,Col_All,Col_temp);
        Pval_All    = cat(1,Pval_All,Pval_temp);
        %disp(['Dist: ' num2str(numel(Dist_All)) ' ,Pval: ' num2str(numel(Pval_All))])
    end
    %% plot V1 effect change
    PDC_All_BR = PDC_All- mean(PDC_All(:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.6),3);
    %
    for t = 1:size(PDC_All_BR,3)
        for f = 1:size(PDC_All_BR,2)
            %[R_all(t,f) P_all(t,f)] = corr(PDC_All_BR(Pval_All>0,f,t),Dist_All(Pval_All>0));
            [R_diagrm(t,f) P_diagrm(t,f)] = corr(PDC_All_BR((Diag_All==0) & (Pval_All>0) & (Col_All==1),f,t),Dist_All((Diag_All==0) & (Pval_All>0)& (Col_All==1)));
        end
    end
    %
    FIG = figure;
    set(FIG,'color','white','unit','inch','position',[5 1 12 4])
    IM = imagesc(StokALL.S766640955.Times,[],R_diagrm');
    set(IM, 'alphadata',P_diagrm'<.01)
    axis xy
    xlim([0 2])
    colormap('jet')
    caxis([-.2 .2])
    colorbar;
    title('FB PDC correlation with probe distance')
    export_fig(FIG,fullfile(Path,[CondName '_FB_DistCorr']),'-pdf','-r200');
    
    %% plot all effect change
    PDC_All_BR = PDC_All- mean(PDC_All(:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.6),3);
    %
    for t = 1:size(PDC_All_BR,3)
        for f = 1:size(PDC_All_BR,2)
            %[R_all(t,f) P_all(t,f)] = corr(PDC_All_BR(Pval_All>0,f,t),Dist_All(Pval_All>0));
            [R_diagrm(t,f) P_diagrm(t,f)] = corr(PDC_All_BR((Diag_All==0) & (Pval_All>0),f,t),Dist_All((Diag_All==0) & (Pval_All>0)));
        end
    end
    %
    FIG = figure;
    set(FIG,'color','white','unit','inch','position',[5 1 12 4])
    IM = imagesc(StokALL.S766640955.Times,[],R_diagrm');
    set(IM, 'alphadata',P_diagrm'<.01)
    axis xy
    xlim([0 2])
    colormap('jet')
    caxis([-.2 .2])
    colorbar;
    title('Golbal PDC correlation with probe distance')
    export_fig(FIG,fullfile(Path,[CondName '_Global_DistCorr']),'-pdf','-r200');
    %% plot between-area effect change
    PDC_All_BR = PDC_All- mean(PDC_All(:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.6),3);
    %
    for t = 1:size(PDC_All_BR,3)
        for f = 1:size(PDC_All_BR,2)
            %[R_all(t,f) P_all(t,f)] = corr(PDC_All_BR(Pval_All>0,f,t),Dist_All(Pval_All>0));
            [R_diagrm(t,f) P_diagrm(t,f)] = corr(PDC_All_BR((Diag_All==0) & (Pval_All>0) & (Laminar_All==0),f,t),Dist_All((Diag_All==0) & (Pval_All>0)& (Laminar_All==0)));
        end
    end
    %
    FIG = figure;
    set(FIG,'color','white','unit','inch','position',[5 5 12 4])
    IM = imagesc(StokALL.S766640955.Times,[],R_diagrm');
    set(IM, 'alphadata',P_diagrm'<.01)
    axis xy
    xlim([0 2])
    colormap('jet')
    %caxis([-.2 .2])
    colorbar;
    title('inter-area PDC correlation with probe distance')
    export_fig(FIG,fullfile(Path,[CondName '_InterArea_DistCorr']),'-pdf','-r200');
    %% plot laminar effect
    PDC_All_BR = PDC_All- mean(PDC_All(:,:,StokALL.S766640955.Times<0 & StokALL.S766640955.Times>-.6),3);
    %
    for t = 1:size(PDC_All_BR,3)
        for f = 1:size(PDC_All_BR,2)
            %[R_all(t,f) P_all(t,f)] = corr(PDC_All_BR(Pval_All>0,f,t),Dist_All(Pval_All>0));
            [R_diagrm(t,f) P_diagrm(t,f)] = corr(PDC_All_BR((Diag_All==0) & (Pval_All>0)& (Laminar_All==1),f,t),Dist_All((Diag_All==0) & (Pval_All>0)& (Laminar_All==1)));
        end
    end
    %
    FIG = figure;
    set(FIG,'color','white','unit','inch','position',[5 11 12 4])
    IM = imagesc(StokALL.S766640955.Times,[],R_diagrm');
    set(IM, 'alphadata',P_diagrm'<.01)
    axis xy
    xlim([0 2])
    colormap('jet')
    title('Laminar PDC correlation with probe distance')
    %caxis([-.5 -.05])
    colorbar;
    export_fig(FIG,fullfile(Path,[CondName '_Laminar_DistCorr']),'-pdf','-r200');
end