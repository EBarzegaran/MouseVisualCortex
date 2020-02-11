function Correlate_hierarchy_plot(Data,Time,Colors)

% Data is a matrix of size: sample x ROI
load Hierarchyscores;
HS = mean(H(Time>-.3 & Time<0,:));

nROIs = size(Data,2);

Ys = [];
Xs = [];
for roi = 1:nROIs
    scatter(repmat(HS(roi),[1 11]), Data(:,roi),50,Colors(roi,:),'filled');hold on;% color and filled, fit line, correlation
    Xs = [Xs squeeze(repmat(HS(roi),[1 11]))];
    Ys = [Ys squeeze(Data(:,roi))'];
end
xlim([-.2 .22])
Xs(isnan(Ys))=[];
Ys(isnan(Ys))=[];
[r,p] = corr(Xs',Ys');
coefficients = polyfit(Xs', Ys', 1);
yFit = polyval(coefficients , -.18:.01:.2);
hold on;
plot(-.18:.01:.2,yFit,'k--','linewidth',1.5);
xlabel('Hierarchy Score');
set(gca,'fontsize',18)
title(['r = ' num2str(round(r,2)) ', p = ' num2str(round(p,5))])
end