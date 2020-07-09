function plot_graph(A,Nodelabels,Colors,direction,Percent,Coords)

% A is a node(from/source) X node(to/target) matrix
A = A';
NNode = size(A,1);
if ~exist('Coords','var')
    if strcmp(direction,'down')
        Coords = NNode:-1:1;
    else
        Coords = (1:NNode);
    end
end


% Node sizes
SNode = A(1:length(A)+1:end);% diagonals
SNode = SNode./max(SNode);
SNode(isnan(SNode))=1;
A(1:length(A)+1:end) = 0;

hold on;
TH = quantile(A(:),Percent);

%disp(sum(A(:)>TH))
%% Plot the edges
NE = 0;
for N1 = 1:NNode
    for N2 = NNode:-1:N1+1
        if A(N1,N2)>TH 
            Y = linspace(Coords(N1),Coords(N2),50);
            X1 = [0, (Coords(N1)-Coords(N2))/4 0];
            X = spline([Coords(N1) (Coords(N1)+Coords(N2))/2 Coords(N2)],X1,Y);
            plot(X,Y,'linewidth',A(N1,N2)*8,'color',Colors(N1,:))
            NE = NE+1;
        end
    end
end
for N2 = NNode:-1:1
    for N1 = 1:N2-1
        if A(N2,N1)>TH 
            Y = linspace(Coords(N2),Coords(N1),50);
            X1 = [0, (Coords(N2)-Coords(N1))/4 0];
            X = spline([Coords(N2) (Coords(N1)+Coords(N2))/2 Coords(N1)],X1,Y);
            plot(X,Y,'linewidth',A(N2,N1)*5,'color',Colors(N2,:))
            NE = NE+1;
        end
    end
end

%% plot the nodes
XL = -min(Coords)+max(Coords);

for N = 1:NNode
    scatter(0, Coords(N),200*SNode(N),Colors(N,:),'filled','MarkerEdgeColor','k');
end

for N = 1:NNode
    if ~isempty(Nodelabels)
        %text(0,Coords(N)+.05*XL,Nodelabels{N},'HorizontalAlignment','center','fontsize',11,'color','k','fontweight','bold');
        text(-max(Coords)/4,Coords(N),Nodelabels{N},'HorizontalAlignment','center','fontsize',14,'color','k','fontweight','normal');
    end
end

%ylim([-1 NNode]+1)
ylim([min(Coords) max(Coords)]*1);
XL = -min(Coords)+max(Coords);
xlim([-XL/2 XL/2]);
set(gca,'xtick',[],'ytick',sort(Coords,'ascend'));%,'yticklabel',[]);
box on
%disp(NE)

end