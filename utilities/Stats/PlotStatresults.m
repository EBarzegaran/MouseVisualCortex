function FIG = PlotStatresults(P,Tval,avdims, time, freq,varargin)


% INPUTs:
    % P
    % Tval
    % avdims
    % time
    % freq
    % <OPTIONS>

%% Default values
opt = ParseArgs(varargin,...
    'figpath'       ,[],...
    'figtitle'      ,'untitiled',...
    'PThresh'       ,.01,... % threshold on P-value
    'SThresh'       ,15,... % color threshold on the test statistics %%%% UPDATE BASED ON THE DATA
    'Twin'          ,[-.3 1],...
    'Labels'        ,arrayfun(@(x) ['L' num2str(x)],1:10,'uni',false),...
    'ColorM'        ,jmaColors('coolhotcortex'),...
    'newfig'        ,true,...
    'colorb'        ,true,...
    'xylabels'      ,true...
    );

%%
FS = 12;

Xticklabels = -.2:.200:1.000;
Xticks = arrayfun(@(x) find(round(time,2)==x,1),Xticklabels);
Yticklabels = 10:20:100;
Yticks = arrayfun(@(x) find(round(freq,2)==x,1),Yticklabels);

pd = setdiff([1 2],avdims);
if opt.newfig
    FIG = figure;
else
    FIG=[];
end
%-----------------------PLOT one-dimension AVERAGED RESULTS---------------    
if numel(pd)==1 % averaged over a dimension
    for i = 1:size(Tval,pd)
        subplot(size(Tval,pd),1,i)
        if pd==1
            %h = imagesc(squeeze(mean(ci(i,:,:,:,:),5)));
            h = imagesc(squeeze(Tval(i,:,:,:)));
        else
            %h = imagesc(squeeze(mean(ci(:,i,:,:,:),5)));
            h = imagesc(squeeze(Tval(:,i,:,:)));
        end
        caxis([-opt.SThresh opt.SThresh])
        colormap(opt.ColorM)
        xlim([find(round(time,2)==opt.Twin(1),1) find(round(time,2)==opt.Twin(2),1)])
        axis xy
        if pd==1
            set(h, 'AlphaData', squeeze(P(i,:,:,:,:)<opt.PThresh)*.9+.1);
        else
            set(h, 'AlphaData', squeeze(P(:,i,:,:,:)<opt.PThresh)*.9+.1);
        end
        title(opt.Labels{i});
        set(gca,'xtick',Xticks,'xticklabel',Xticklabels*1000,'ytick',Yticks,'yticklabel',Yticklabels,'fontsize',16);
        if i==size(Tval,pd)
            ylabel('Frequency (Hz)');
            xlabel('Time (msec)');
            ps = get(gca,'position');
%             cb = colorbar;
%             set(gca,'position',ps);
%             set(cb,'position',get(cb,'position')+[-.05 0 0 0])
%             set(get(cb,'title'),'string','T-Statistics');
        end
        colorbar

    end
    dims = size(P);
    dims = sort(dims(1:2));
    if opt.newfig
        set(FIG,'unit','inch','position',[1 1 dims(1)*9 dims(2)*3.5],'color','w');
    end
%-----------------------PLOT AVERAGED RESULTS-----------------------------        
elseif numel(pd)==0
    h = imagesc(time(1:end-5),[],squeeze(Tval(:,:,:,1:end-5,:)));
    caxis([-opt.SThresh opt.SThresh])
    %caxis([0 opt.SThresh])
    colormap(opt.ColorM)
    %xlim([find(round(time,2)==opt.Twin(1),1) find(round(time,2)==opt.Twin(2),1)])
    xlim([opt.Twin(1) opt.Twin(2)])
    axis xy
    set(h, 'AlphaData', squeeze(P(:,:,:,1:end-5,:)<opt.PThresh)*.9+.1);
    title(opt.figtitle);
    set(gca,'xtick',Xticklabels,'xticklabel',Xticklabels,'ytick',Yticks,'yticklabel',Yticklabels,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.025]);
    box off;
    if opt.xylabels
        xlabel('Time(S)');
        ylabel('Frequency(Hz)');
    end
    if opt.colorb
        
        GCP = get(gca,'position');
        CB = colorbar;
        colorTitleHandle = get(CB,'Title');
        set(colorTitleHandle ,'String','T-Value');
        set(gca,'position',GCP);
        set(CB, 'ylim', [0 opt.SThresh])
    end
    set(gca,'fontsize',FS)
    dims = size(P);
    dims = sort(dims(1:2));
    if opt.newfig
        set(FIG,'unit','inch','position',[1 1 dims(1)*8 dims(2)*3.5],'color','w');
    end
%-----------------------PLOT NON AVERAGED RESULTS-------------------------       
elseif numel(pd)==2
    for i = 1:size(Tval,1)
        for j = 1:size(Tval,2)
            subplot(size(Tval,1),size(Tval,2),(i-1)*size(Tval,1)+j);
            h = imagesc(squeeze(Tval(i,j,:,:,:)));
            caxis([-opt.SThresh opt.SThresh])
            colormap((opt.ColorM));
            xlim([find(round(time,2)==opt.Twin(1),1) find(round(time,2)==opt.Twin(2),1)])
            axis xy
            set(h, 'AlphaData', squeeze(P(i,j,:,:,:)<opt.PThresh)*.8+.2);
            %----------axis adjust------------------
            if i==1, title(opt.Labels{j});end
            if j==1 
                ylabel(['L' num2str(i)],'fontweight','bold');
                set(gca,'ytick',Yticks(1:2:end),'yticklabel',Yticklabels(1:2:end));
            else
                set(gca,'ytick',Yticks,'yticklabel',[]); 
            end
            if i==size(Tval,1)
                
                if j==size(Tval,2)
                    set(gca,'xtick',Xticks(1:2:end),'xticklabel',Xticklabels(1:2:end)*1000);
                    xlabel('Time (S)');
                    ylabel('Freq (Hz)');
                    ps = get(gca,'position');
                    colorbar;
                    set(gca,'position',ps);
                else
                    set(gca,'xtick',Xticks,'xticklabel',[]);
                end
            else
                set(gca,'xtick',Xticks,'xticklabel',[]);
            end
            set(gca,'fontsize',16);
        end
    end
    dims = size(P);
    dims = sort(dims(1:2));
    if opt.newfig
        set(FIG,'unit','inch','position',[1 1 dims(1)*5 dims(2)*3.5],'color','w');
    end
end
%--------------figure title-------------------------
%axes('position',[.5 .98 .1 .05]); axis off;
%text(0,0,opt.figtitle,'HorizontalAlignment','center','fontsize',16);
if opt.newfig
    export_fig(fullfile(opt.figpath,opt.figtitle),'-pdf','-r200');close;
end
end