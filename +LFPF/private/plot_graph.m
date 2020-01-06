function plot_graph(A,Nodelabels,Colors,direction,Percent)

% A is a node(from/source) X node(to/target) matrix
A = A';
NNode = size(A,1);
if strcmp(direction,'down')
    Coords = NNode:-1:1;
else
    Coords = (1:NNode);
end
hold on;
TH = quantile(A(:),Percent);

%% Plot the edges
for N1 = 1:NNode
    for N2 = NNode:-1:N1+1
        if A(N1,N2)>TH 
            Y = linspace(Coords(N1),Coords(N2),50);
            X1 = [0, (Coords(N1)-Coords(N2))/4 0];
            X = spline([Coords(N1) (Coords(N1)+Coords(N2))/2 Coords(N2)],X1,Y);
            plot(X,Y,'linewidth',A(N1,N2)*5,'color',Colors(N1,:))
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
        end
    end
end

%% plot the nodes
for N = 1:NNode
    scatter(0, Coords(N),50,Colors(N,:),'filled');
    if ~isempty(Nodelabels)
        text(0,Coords(N)+.2,Nodelabels{N},'HorizontalAlignment','center','fontsize',8);
    end
end
ylim([-1 NNode]+1)
xlim([-NNode/3 NNode/3]);
set(gca,'xtick',[],'ytick',[]);
box on

end