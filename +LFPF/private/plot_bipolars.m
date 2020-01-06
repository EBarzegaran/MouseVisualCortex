function plot_bipolars(Data,LayerInfo,Probe_ID,Save_Path)
    
Time = Data.Times;
Colors = brewermap(6,'Spectral');
FIG = figure;
set(FIG,'unit','inch','position',[1 1 8 15],'color','w')
for i = 1:4
    subplot(4,1,i);
    data = squeeze(mean(Data.Y(:,:,:,i),1));
    imagesc(data);
    
    set(gca,'xtick',round(1:numel(Time)/20:numel(Time)),'xticklabels',round(Time(round(1:numel(Time)/20:numel(Time))),2),...
        'ytick',1:4:22,'yticklabel',arrayfun(@(x) num2str(x),(1:4:22)*40,'uni',false));
    vline(min(find(round(Data.Times,2)==0)),'w--');
    xlim([min(find(round(Time,2)==-.1)) min(find(round(Time,2)==.5))]);
    caxis([-max(abs(data(:))) max(abs(data(:)))])
    colorbar;
    title([ 'Orientation = ' num2str(Data.cnd_info.orientation(i))]);
    
    for l = 1:6
        hline(LayerInfo(l),'k--',['L' num2str(l)]);
    end
end
colormap(jmaColors('coolhotcortex'));

export_fig(fullfile(Save_Path,[Probe_ID '_LFP_Bipolar_Gratings']))
close;
end