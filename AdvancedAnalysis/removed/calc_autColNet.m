% visualizing inter-areal and laminar conenctivity in the Allen data
% 
% calculate iPDC with autoregressive parts for all animals
% comparing contributions of the auto-regressive, the within-collum and
% inter-areal coeficients.


%%
%clean cleaner cleanest
clear; clc
%%
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/FullModel/AR_ALL_drifting_gratings_75_repeats__contrast0-8_iPDC_Mord8_ff098.mat');
fields = fieldnames(AR_ALL); % animal ids


%% loop through struct, calc iPDC, separate influences
% aut - autoregressive parts
% net - inflows from entire network
% col - inflows from within the roi, collumn
% ext - inflows from other areas

srate = 250;
freqs = AR_ALL.(fields{1}).f; % start:end freq, Hz
univ = 1;
flow = 1; 

for f = 1:numel(fields)
    disp([num2str(f), ' out of ' num2str(numel(fields))])
    labels = AR_ALL.(fields{f}).ROIs; % area names
    SK = AR_ALL.(fields{f}).KF;
    PDC = dynet_ar2pdc(SK,srate,freqs,'iPDC',univ,flow,0);
    
    dims = size(PDC);
    aut = zeros(dims(2:4)); %autoregressive part
    for dg = 1:dims(1) % extract and clear diagonals
         aut(dg,:,:) = PDC(dg, dg, :, :);
         PDC(dg, dg, :, :) = 0; 
    end

    net = squeeze(sum(PDC, 2)); %all influence from elsewehere, excluding auto

    col = zeros(dims(2:4)); % from within the ROI
    ext = zeros(dims(2:4)); % from other ROIs

    for roi = 0:(numel(labels)-1)
        sIdx = (roi * 6) + 1;
        eIdx = roi * 6 + 6;

        extIdx = 1:dims(1); % all indices
        extIdx(sIdx:eIdx) = []; % idx of elements outside the collumn

        for i = sIdx:eIdx % all layers of the ROI
            col(i,:,:) = sum(PDC(i, sIdx:eIdx, :, :),2); % inflow from colum  
            ext(i,:,:) = sum(PDC(i, extIdx, :, :),2); % inflow from outside the roi      
        end
    end
    
    % store
    AR_ALL.(fields{f}).iPDC = PDC;
    AR_ALL.(fields{f}).aut = aut;
    AR_ALL.(fields{f}).net = net;
    AR_ALL.(fields{f}).col = col;
    AR_ALL.(fields{f}).ext = ext;    
end


%
save('AR_ALL_drifting_gratings_75_repeats__contrast0-1_iPDC_Mord8_ff098.mat', 'AR_ALL', '-V7.3')

%% average across animals
load('AR_ALL_drifting_gratings_75_repeats__contrast0-1_iPDC_Mord8_ff098.mat');

ROIs_Select ={'VISp','VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'}; % ordered
% Path = '';
% SaveName = 'iPDC_ff.98_MOrd15_17022020.mat';
% STOKAllAverage(AR_ALL, ROIs_Select, Path,SaveName)


%% plot all
% aut - autoregressive parts
% net - inflows from entire network
% col - inflows from within the roi, collumn
% ext - inflows from other areas

f=8;
labels = AR_ALL.(fields{f}).ROIs; % area names
msec=round((AR_ALL.(fields{f}).Times)*1000); 
freqs = AR_ALL.(fields{1}).f; % start:end freq, Hz
aut = AR_ALL.(fields{f}).aut;
net = AR_ALL.(fields{f}).net;
col = AR_ALL.(fields{f}).col;
ext = AR_ALL.(fields{f}).ext;


nAreas = numel(labels);
close all % plot

% limits
stf = find(msec > -50, 1);
etf= find(msec > 1900, 1);
sfp = 1 ;
efp = 100;


%dat= aut-col;

% proportional
%dat= (aut-net)./aut; % aut - net
%dat= (aut-col)./aut; % aut - col
%dat= (aut-ext)./aut; % aut - ext
%dat= (col-ext)./col; % col - ext


% percentage of total
dat1 = (aut);%./(ext+col+aut)) *100;
dat2 = (col);%./(ext+col+aut)) *100;
dat3 = (ext);%./(ext+col+aut)) *100;


%zmax = max(abs(max(dat(:))), abs(min(dat(:)))); % 
%zmin = -zmax; 
zmin = 0; % percentage
zmax = 100;

figure('units','normalized','position',[0 0 .9 .9]);
colormap(jet)
%colormap(hot) % percentage

fig=1;
for l = 1:6 % layers    
    for r = 1:nAreas % areas
        idx = (r * 6) - 6 + l; % area, layer
        subplot(6,nAreas,fig)
        hold on;
        
        datN = dat1(idx,:,stf:etf); %- nanmean(dat(idx,:,msec>-300 & msec<0),3);
        plot( freqs, squeeze(mean(datN,3)),'linewidth',1.5);%-aut(idx,:,:)+ext(idx,:,:)));
        
        datN = dat2(idx,:,stf:etf); %- nanmean(dat(idx,:,msec>-300 & msec<0),3);
        plot( freqs, squeeze(mean(datN,3)),'linewidth',1.5);%-aut(idx,:,:)+ext(idx,:,:)));
       
        datN = dat3(idx,:,stf:etf); %- nanmean(dat(idx,:,msec>-300 & msec<0),3);
        plot( freqs, squeeze(mean(datN,3)),'linewidth',1.5);%-aut(idx,:,:)+ext(idx,:,:)));

        %imagesc(msec, freqs, squeeze(ext(idx,:,:)+aut(idx,:,:)+col(idx,:,:)));
        axis xy
        vline(0)
        title(strcat(string(labels(r)), ' L', num2str(l)))
        if fig==1
            legend({'aut','col','ext'});
        end
        
        fig=fig+1;
        
    end
end

colorbar('Location','East')
dynet_SSM_KF

%%
% conclusion: The autoregressive part is looking weird!



