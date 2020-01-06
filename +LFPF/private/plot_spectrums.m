function FIG = plot_spectrums(RS2,Times,F, Nodes,ElecLabels,Title)
% This function plots the time-frequency spectrum of the data, with their
% marginal distribution

%%
MaxS = max(max(max(RS2(:,15:end,:))));
%
Nodes = min(Nodes,6);

num_time_labels = 6;
time_label_indices = round(linspace(1, numel(Times), num_time_labels));
time_labels = round(Times(time_label_indices),2);

num_F_labels = 7;
F_label_indices = round(linspace(1, numel(F), num_F_labels));
F_labels = round(F(F_label_indices),2);


FIG = figure;
set(FIG,'unit','inch','position',[0 0 25 10],'color','w')
for E = 1:Nodes
    SP(E) = subplot(2,3,E); hold on;
    imagesc(abs(squeeze(RS2(E,:,:))));
    %---------------------x and y ticks--------------------------
    if E ==7
        xlabel('Time(s)');        
    end
    set(gca,'ytick',F_label_indices,'yticklabel',[],'xtick',time_label_indices,'xticklabel',time_labels);

    %---------------------lines----------------------------------
    for yl = F_label_indices
        line([0 numel(Times)],[yl yl],'linestyle','--','color','w')
    end
    TimeZero = find(round(Times,3)==0);
    line([TimeZero(1) TimeZero(1)],[1 100],'linestyle','--','color','w','linewidth',2);
    for xl = 1:numel(time_label_indices)
        vline(time_label_indices(xl),'w--');
    end
    %--------------axis control----------------------------------
    axis xy tight;
    %caxis([0 MaxS/2])
    ylim([1 numel(F)]);
    TimeInit = find(round(Times,3)==-.1);
    %xlim([TimeInit(1) min(find(round(Times,2)==.2))]);%numel(Times)]);
end


colormap('jet')
%-----------------adjust subplots and add marginal distributions-----------
Offvalues = mod(1:Nodes,3);
Offvalues(Offvalues==0)=3;
for E = 1:Nodes
    set(SP(E),'position',get(SP(E),'position')-[(-.019*Offvalues(E))+.07 0.0 0 .05])
    Pos = get(SP(E),'position');
    %-----------------------------spectrum---------------------------------
    A = axes('position',[Pos(1)-.03 Pos(2) .03 Pos(4)] ); %axis off
    plot(A,squeeze(abs(mean(RS2(E,:,Times>0),3))),'linewidth',2,'color','k')
    view([-90 90])
    xlim([1 numel(F)])
    if mod(E,3)==1
        set(gca,'xtick',F_label_indices,'xticklabel',F_labels,'ytick',[]);
    else
        set(gca,'xtick',F_label_indices,'xticklabel',[],'ytick',[]);
    end
    
    for yl = F_label_indices
        vline(yl,'k--')
    end
    if E==7
        xlabel('Frequerncy(Hz)');
    end
    %---------------------------TimeAverage--------------------------------
    A = axes('position',[Pos(1) Pos(2)+Pos(4) Pos(3) .07] ); %axis off
    hold on;
    %
    plot(A,squeeze(mean(abs(RS2(E,F<=20,:)),2)),'linewidth',1.5,'color','b');
    plot(A,squeeze(mean(abs(RS2(E,F>20 & F<=40,:)),2)),'linewidth',1.5,'color','r');
    plot(A,squeeze(mean(abs(RS2(E,F>40,:)),2)),'linewidth',1.5,'color','k');
    set(gca,'xtick',F_label_indices,'xticklabel',[],'ytick',[]); box off    
    for xl = 1:numel(time_label_indices)
        vline(time_label_indices(xl),'k--');
    end
    if E ==Nodes
        %L = legend('[1,80] Hz','[5,30] Hz','[30,80] Hz');
        L = legend('<20 Hz','>20 & <40 Hz','>40Hz');
        set(L,'position',get(L,'position')+[0.02 0 0 0]);
    end
   xlim([1 numel(Times)]);%numel(Times)]);
    %--------------------title-----------------------------------
    title([ElecLabels{E} ' (' num2str((E-1) * 40 *2) '\mu m)'])
end

%% Title

if exist('Title','var') & ~isempty(Title)
    T_AX = axes('position',[.4 .95 .2 .05]);
    axis off;
    text(0.5,.8,Title,'HorizontalAlignment','center','fontsize',12,'fontweight','bold')
end
end