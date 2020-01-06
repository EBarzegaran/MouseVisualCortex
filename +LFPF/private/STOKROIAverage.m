function [STOK_Averaged, FIG] = STOKROIAverage(STOK_Res,ROIselect,savefig,Path,PDCMethod)

Sessions_ID = fieldnames(STOK_Res);

ROIs = cellfun(@(x) fieldnames(STOK_Res.(x)),Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(1,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

for roi = 1:numel(ROIs)
    for S = 1:numel(Sessions_ID)
        if isfield(STOK_Res.(Sessions_ID{S}),ROIs{roi})
            Temp{S} =  STOK_Res.(Sessions_ID{S}).(ROIs{roi}).PDC;
            STOK_Averaged.(ROIs{roi}).Times = STOK_Res.(Sessions_ID{S}).(ROIs{roi}).Times;
            STOK_Averaged.(ROIs{roi}).f = STOK_Res.(Sessions_ID{S}).(ROIs{roi}).f;
        end
    end
    STOK_Averaged.(ROIs{roi}).PDC = mean(cat(5,Temp{:}),5);
    STOK_Averaged.(ROIs{roi}).Num = ROI_counts(roi);
    
    %---------------------plot the spectrum--------------------------------
    Nodes = size(STOK_Averaged.(ROIs{roi}).PDC,1);
    TimesT = STOK_Averaged.(ROIs{roi}).Times;
    
    Temp_PSD = STOK_Averaged.(ROIs{roi}).PDC;
    Temp_PSD = reshape(Temp_PSD,Nodes*Nodes,size(Temp_PSD,3),size(Temp_PSD,4));
    Temp_PSD =  Temp_PSD(1:Nodes+1:end,:,:);
    Temp_PSD =  Temp_PSD - mean(Temp_PSD(:,:,TimesT<0 & TimesT>-.2),3);
    
    Title = [ROIs{roi} ' - Averaged over ' num2str(STOK_Averaged.(ROIs{roi}).Num) ' sessions'];
    FIG = plot_spectrums(Temp_PSD(:,STOK_Averaged.(ROIs{roi}).f,TimesT>-.1 & TimesT<1),TimesT(TimesT>-.1 & TimesT<1), STOK_Averaged.(ROIs{roi}).f(10:80),Nodes,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false),Title);
    export_fig(FIG,fullfile(Path,['PSD_individual_' ROIs{roi}]),'-pdf');
    close;
end
clear FIG;
%% plot the results
if ~isempty(ROIselect)
    ROIs = intersect(ROIs,ROIselect);
end

for roi = 2:numel(ROIs)
    S_STOK = STOK_Averaged.(ROIs{roi}).PDC;
    Times = STOK_Averaged.(ROIs{roi}).Times;
    f = STOK_Averaged.(ROIs{roi}).f;
    S_STOKN     = S_STOK - mean(S_STOK(:,:,:,(Times<0)& (Times>-0.3)),4);
    Title = [ROIs{roi} ' - Averaged over ' num2str(STOK_Averaged.(ROIs{roi}).Num) ' sessions'];
    FIG{roi} = dynet_connplot(S_STOKN(:,:,:,Times>-.2 & Times<.3) ,Times(Times>-.1 &  Times<.2),f,[],[.05 .95],[],[],1,Title);
    set(FIG{roi},'color','w','unit','inch','position',[0 0 25 15])
    if savefig
        export_fig(FIG{roi},fullfile(Path,['STOK_individual_Transient_' ROIs{roi} '_' PDCMethod]),'-pdf');
    end
end

end