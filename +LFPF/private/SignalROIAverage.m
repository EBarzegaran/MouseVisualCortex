function [Signal_Averaged, FIG] = SignalROIAverage(Session,ROIselect,savefig,Path,SaveName,SpecEstim,Freqs)

Sessions_ID = fieldnames(Session);

ROIs = cellfun(@(x) Session.(x).ROIs,Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(2,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

for roi = 1:numel(ROIs)
    disp(ROIs{roi})
    for S = 1:numel(Sessions_ID)
            
        if isfield(Session.(Sessions_ID{S}),ROIs{roi})
            t = Session.(Sessions_ID{S}).(ROIs{roi}).Times;
            if (S==1) && (roi==1)
                Ind = 1:numel(t);
            end
            
            Temp{S} =  squeeze(mean(mean(Session.(Sessions_ID{S}).(ROIs{roi}).Y(:,:,Ind,:),1),4));
            P(S)    =  mean((mean(Temp{S}(:,t>0).^2,2)),1)*250;
            SNR(S)  =  max((max(Temp{S}(:,t>0).^2,[],2)),[],1)./mean((mean(Temp{S}(:,t<0 & t>-.3).^2,2)),1);
            
            % FFT
            Session.(Sessions_ID{S}).(ROIs{roi}).Y = Session.(Sessions_ID{S}).(ROIs{roi}).Y(:,:,Ind,:);
            Session.(Sessions_ID{S}).(ROIs{roi}).Times = Session.(Sessions_ID{S}).(ROIs{roi}).Times(Ind);
            [PostStim,PreStim]  =  FFT_estimate(Session.(Sessions_ID{S}).(ROIs{roi}),1);
            
            FFT{S} = PostStim.FFT;
            AMP{S} = PostStim.AMP;
            F      = PostStim.F;
            
            FFT_Pres{S} = PreStim.FFT;
            AMP_Pres{S} = PreStim.AMP;
            F_Pres      = PreStim.F;
            if S <5
                Signal_Averaged.(ROIs{roi}).Times = Session.(Sessions_ID{S}).(ROIs{roi}).Times(Ind);
            end
        end
    end
    % Signal Data
    Signal_Averaged.(ROIs{roi}).SM  = mean(cat(3,Temp{:}),3);
    Signal_Averaged.(ROIs{roi}).STD = std(cat(3,Temp{:}),[],3);
    Signal_Averaged.(ROIs{roi}).Num = ROI_counts(roi);
    Signal_Averaged.(ROIs{roi}).P   = P;
    Signal_Averaged.(ROIs{roi}).SNR   = SNR;
    % FFT Data
    Signal_FFT_Averaged.(ROIs{roi}).SMEvoked  = cat(3,FFT{:});
    Signal_FFT_Averaged.(ROIs{roi}).SMSpon  = cat(3,AMP{:});
    Signal_FFT_Averaged.(ROIs{roi}).Num = ROI_counts(roi);
    Signal_FFT_Averaged.(ROIs{roi}).Freq = F;
    Signal_FFT_Averaged.(ROIs{roi}).FreqPres = F_Pres;
    Signal_FFT_Averaged.(ROIs{roi}).SMPres = cat(3,AMP_Pres{:});
    Signal_FFT_Averaged.(ROIs{roi}).SMPresPhase = cat(3,FFT_Pres{:});
    clear Temp P SNR FFT;
end
save(fullfile(Path,'Fullmodel',['Signal_Averaged' SaveName]),'Signal_Averaged','Signal_FFT_Averaged','-v7.3');

%% Spectrum estimation
if SpecEstim
    if ~exist(fullfile(Path,'Fullmodel',['Signal_PSD' SaveName '.mat']))
        for roi = 1:numel(ROIs)
            disp(ROIs{roi})
            for S = 1:numel(Sessions_ID)
                S
                if isfield(Session.(Sessions_ID{S}),ROIs{roi})
                    %[Wavelet_Temp{S} STOK_Temp{S}] = LFPF.SpecEstimate(Session.(Sessions_ID{S}).(ROIs{roi}),Freqs);
                    Data = Session.(Sessions_ID{S}).(ROIs{roi});
                    Data.Y = Data.Y(:,:,Ind,:);
                    [Wavelet_Temp{S} STOK_Temp{S}] = LFPF.SpecEstimate(Data,Freqs);
                else
                    Wavelet_Temp{S} =[];
                    STOK_Temp{S}    =[];
                end
            end
            %WaveletPSD.(ROIs{roi}) = mean((cat(4,Wavelet_Temp{:})),4);
            WaveletPSD.(ROIs{roi}) = Wavelet_Temp;
            STOKPSD.(ROIs{roi}) = STOK_Temp;
            clear Wavelet_Temp STOK_Temp;
        end
        save(fullfile(Path,'Fullmodel',['Signal_PSD' SaveName]),'WaveletPSD','STOKPSD','t','ROIs','Freqs','-v7.3');
    else
        load(fullfile(Path,'Fullmodel',['Signal_PSD' SaveName '.mat']));
    end
end

%% plot the time domain results
% first the order of ROIs

Colors = brewermap(6,'Spectral');
if ~isempty(ROIselect)
    ROIs = intersect(ROIs,ROIselect);
    ord = cellfun(@(x) find(strcmp(ROIs,x)),ROIselect);
    ROIs = ROIs(ord);
end

% load ROI names for figures
load ROInames;

COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs,'SubColors');
%---------------------------------------------------
Times = Signal_Averaged.(ROIs{roi}).Times*1000;
Xticks = arrayfun(@(x) find(round(Signal_Averaged.(ROIs{roi}).Times,2)*1000==x,1),[ 0:200:1000])+1;

Yticks = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
%--------------------------------------------------

FIG = figure;

for roi = 1:numel(ROIs)
    
    
    
    if numel(ROIs)>6
        subplot(ceil(numel(ROIs)/2),2,roi);
    else
        subplot(numel(ROIs),1,roi)
    end
    hold on;
    Signal = Signal_Averaged.(ROIs{roi}).SM;
    SignalErr = Signal_Averaged.(ROIs{roi}).STD;
    if roi==1
        Offset = max(Signal(:))*1;
    end
    for l = 1:size(Signal,1)
        X = Signal_Averaged.(ROIs{roi}).Times*1000;
        Y = Signal(l,:)-Offset*l;
        Err = SignalErr(l,:)/sqrt(Signal_Averaged.(ROIs{roi}).Num);
       % fill([X X(end:-1:1)],[Y+Err Y(end:-1:1)-Err(end:-1:1)],Colors(roi,l,:),'facealpha',.3,'edgecolor','none');
    end
    
    for l = 1:size(Signal,1)
        LG(l)=plot(Signal_Averaged.(ROIs{roi}).Times*1000,Signal(l,:)-Offset*l,'linewidth',1.5,'color',Colors(roi,l,:));
        MS(l) = mean(Signal(l,:)-Offset*l);
    end
    xlim([-100 1001])
    %xlim([40 401])
    
    %title([ROI_names.(ROIs{roi}) ' - Averaged over ' num2str(Signal_Averaged.(ROIs{roi}).Num) ' sessions'])
    TL = title(ROI_names.(ROIs{roi}));
    set(TL,'position',get(TL,'position')+[0 -.0001 0])
    
    if roi == numel(ROIs) 
        xlabel('Time (S)')
         %set(gca,'ytick',MS(end:-1:1),'yticklabel',Yticks(end:-1:1));
    else
        %set(gca,'ytick',MS(end:-1:1),'yticklabel',[]);
    end
    vline(0.0,'k--')
    grid on;
  
    set(gca,'ytick',MS(end:-1:1),'yticklabel',Yticks(end:-1:1),'xtick',Times(Xticks),'xticklabel',round(Times(Xticks))/1000,'fontsize',14,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035]);
    %set(gca,'position',get(gca,'position')+[0 0 0 0.02])
%     if roi ==2
%         L = legend(LG,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
%         set(L,'position',get(L,'position')+[.1 0 0 0])
%         
%     end
end

set(FIG,'color','w','unit','inch','position',[0 0 5 16])
if savefig
  export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_transient' SaveName]),'-pdf');
end

%% plot the time domain results
Colors = COBJ.MatrixColors(ROIs,'subColors');
FIG = figure;

for roi = 1:numel(ROIs)
    if roi==1
        M = max(abs(Signal_Averaged.(ROIs{roi}).SM(:)));
    end
    subplot(numel(ROIs),1,roi);
    for l = 1:6
        plot(Times,(Signal_Averaged.(ROIs{roi}).SM(l,:)),'Color',squeeze(Colors(roi,l,:)),'linewidth',1.5);
        hold on;
    end
    [CSD(roi)] = sqrt(mean(mean(Signal_Averaged.(ROIs{roi}).SM(:,Times>0 & Times<2000).^2)));
    hline(0,'k--')
    xlim([-100 1001]);
    ylim([-M M])
    grid on;
    title(ROI_names.(ROIs{roi}));
    if roi== numel(ROIs)
        xlabel('Time (msec)');
        ylabel('Amplitude (\muV)')
    end
    set(gca,'fontsize',16,'xtick',Times(Xticks),'xticklabel',round(Times(Xticks)))
end
set(FIG,'color','w','unit','inch','position',[0 0 5 16])
if savefig
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_Overlap2' SaveName]),'-pdf');
end
close all
%legend(ROIs)
%% plot the frequency domain results
FIG1 = figure(1);
set(FIG1,'color','w','unit','inch','position',[0 0 20 15])

FIG2 = figure(2);
set(FIG2,'color','w','unit','inch','position',[0 0 20 15]);

for roi = 1:numel(ROIs)
    FFT_evoked  = squeeze(abs(mean(Signal_FFT_Averaged.(ROIs{roi}).SMEvoked,3)));
    MEvoked     = max(FFT_evoked(:));
    mEvoked     = min(FFT_evoked(:));
    
    FFT_spon  = squeeze(abs(mean(Signal_FFT_Averaged.(ROIs{roi}).SMSpon,3)));
    MSpon     = max(FFT_spon(:));
    mSpon     = min(FFT_spon(:));
    
    for i = 1:6
        figure(1)
        subplot(6,numel(ROIs),(i-1)*numel(ROIs)+roi)
        plot(Signal_FFT_Averaged.(ROIs{roi}).Freq,FFT_evoked(i,:,:),'linewidth',1.5,'color',squeeze(Colors(roi,i,:)))
        xlim([1 80])
        ylim([mEvoked MEvoked])
         if i==1
            title(ROI_names.(ROIs{roi}));
        end
        if roi==1
            ylabel(['L' num2str(i)]);
        end
        if (roi==numel(ROIs)) && (i==6)
            
            xlabel('Frequency (Hz)');
        end
        
        figure(2)
        subplot(6,numel(ROIs),(i-1)*numel(ROIs)+roi)
        plot(Signal_FFT_Averaged.(ROIs{roi}).Freq,FFT_spon(i,:,:),'linewidth',1.5,'color',squeeze(Colors(roi,i,:)))
        xlim([1 80])
        ylim([mSpon MSpon])
        
        if i==1
            title(ROI_names.(ROIs{roi}));
        end
        if roi==1
            ylabel(['L' num2str(i)]);
        end
        if (roi==numel(ROIs)) && (i==6)
            
            xlabel('Frequency (Hz)');
        end
    end
end

if savefig
    export_fig(FIG1,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Evoked' SaveName]),'-pdf','-d200');
    export_fig(FIG2,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Spontaneous' SaveName]),'-pdf','-d200');
end
% close all;
%% Wavelet results
if SpecEstim
%---------------------------------------------------
Times = Signal_Averaged.(ROIs{roi}).Times*1000;
Xticks = arrayfun(@(x) find(round(Signal_Averaged.(ROIs{roi}).Times,2)*1000==x,1),[ 0:200:1000]);

FIG = figure(1);
set(FIG,'color','w','unit','inch','position',[0 0 10 12])

if false
    FIG2 = figure(2);
    set(FIG2,'color','w','unit','inch','position',[0 0 25 20]);
end
%WaveletPSD = cellfun(@(y) cellfun(@(x) x(:,:,1:376),WaveletPSD.(y),'uni',false),fieldnames(WaveletPSD));

WaveletPSD2 = cellfun(@(x) mean((cat(4,WaveletPSD.(x){:})),4),ROIs,'uni',false);

MWav = max(max(mean(abs(WaveletPSD2{1}))));
for roi = 1:numel(ROIs)
    if false
        figure(2);
        WW = WaveletPSD2{roi};
        WW2 = log10(abs(mean(WaveletPSD2{roi}(:,:,t<0 & t>-.3),3)));
        for i = 1:6
            subplot(6,numel(ROIs),(i-1)*numel(ROIs)+roi)
            plot(log10(Freqs),(log10(abs(squeeze(mean((WW(i,:,t>0)),3))))-squeeze(WW2(i,:))),'k','linewidth',1.5);
            %ylim([0 max(WW(:))*.8])
            ylim([-.2 1])
            xlim([0 log10(100)])
            set(gca,'xtick', log10(Freqs([1 5 10 20 30 50 100])),'xticklabel',Freqs([1 5 10 20 30 50 100]))
            if i==1
                title(ROI_names.(ROIs{roi}));
            end
            if roi==1
                ylabel(['L' num2str(i)]);
            end
            if (roi==numel(ROIs)) && (i==6)
                xlabel('Time (Sec)');
                ylabel('Frequency (Hz)');
            end
        end
    end
    

    figure(1)
    subplot(numel(ROIs),2,(roi-1)*2+2);
    BRWV = abs(squeeze(mean(WaveletPSD2{roi}-(mean(WaveletPSD2{roi}(:,:,t>-.3 & t<0),3))))); %baseline removed
    contourf(t(~isnan(t)),log10(Freqs(Freqs>2)),10*log10(BRWV(Freqs>2,~isnan(t))),100,'linecolor','none');
    caxis([-12 -8]*10);
    set(gca,'ytick',log10([5 10 20 30 50 100]),'yticklabel',[5 10 20 30 50 100])
    %xlim([-.2 1]);
    xlim([.02 .4]);
    
    vline(0,'k--');
    if roi ==numel(ROIs)
        xlabel('Time (Sec)');
        ylabel('Frequency (Hz)');
        SPP = get(gca,'position');
        CB = colorbar;
        set(gca,'position',SPP);
        set(CB,'position',get(CB,'position')+[-.04 0 0 0])
    end
    title(ROI_names.(ROIs{roi}));
    set(gca,'fontsize',14)
    
    
    subplot(numel(ROIs),2,(roi-1)*2+1);
    BRWV = abs(squeeze(mean(WaveletPSD2{roi}-(mean(WaveletPSD2{roi}(:,:,t>-.3 & t<0),3))))); %baseline removed
    contourf(t(~isnan(t)),(Freqs(Freqs>2)),(BRWV(Freqs>2,~isnan(t))),100,'linecolor','none');
    if roi ==1
        MMax = max(BRWV(:))/20;
    end
    caxis([0 MMax]);
    %set(gca,'ytick',[5 10 20 30 50 100],'yticklabel',[5 10 20 30 50 100])
    %xlim([-.2 1]);
    xlim([.02 .4]);
    vline(0,'k--');
    if roi ==numel(ROIs)
        xlabel('Time (Sec)');
        ylabel('Frequency (Hz)');
        SPP = get(gca,'position');
        CB = colorbar;
        set(gca,'position',SPP);
        set(CB,'position',get(CB,'position')+[-.04 0 0 0])
    end

        
    title(ROI_names.(ROIs{roi}));
    set(gca,'fontsize',14)

end
colormap('jet')

if savefig
    
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_PSD_transient' SaveName]),'-pdf','-d200');
    %export_fig(FIG2,fullfile(Path,['Signal_Averaged_AllROIs_PSD_layers' SaveName]),'-pdf','-d200');
end
end
%close all

%% PLOT correlation of signal power with hierarchical scores
Times = Signal_Averaged.(ROIs{roi}).Times*1000;
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
load hierarchyscores;
HM = mean(H(Times>-300 & Times<0,:));

FIG = figure;
Xs = [];
Ys = [];
for roi = 1:numel(ROIs)
    
    P = (Signal_Averaged.(ROIs{roi}).SNR);
    P = real(10*log10(P(P~=0)));
    %std(P)./sum()
    scatter(repmat(HM(roi),1,numel(P)),P,50,Colors(roi,:),'filled');
    %scatter(HM(roi),nanmedian(P),50,Colors(roi,:),'filled');
    Xs = [Xs repmat(HM(roi),1,numel(P))];
    Ys = [Ys P];
    hold on;
end

[R,pval]= corr(Xs',Ys');
coefficients = polyfit(Xs',Ys', 1);
yFit = polyval(coefficients ,-.18:.01:.2);
plot(-.18:.01:.2,yFit,'k--','linewidth',1.5)
[r,p] = corr(Xs',Ys');
xlim([-.2 .22])
%ylim([-87 -65])
set(gca,'fontsize',14)
title(['r = ' num2str(round(r,2)) ', p = ' num2str(round(p,5))])
xlabel('Hierarchy Score');
ylabel('SNR (dB)');
set(FIG,'unit','inch','position',[5 5 5 4],'color','w');
set(gca,'fontsize',18)
if savefig
    export_fig(FIG,fullfile(Path,['Signal_SNR_AllROIs_hierarchyscore' SaveName]),'-pdf');
end
end

