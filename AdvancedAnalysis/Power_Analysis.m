
clear; clc;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';

File1 = 'Signal_Averaged_drifting_gratings_75_repeats__contrast0-8'; %'Signal_Averaged_dot_motion__Speed0-0005-------0-001-------0-005.mat';% 
File2 = 'Signal_Averaged_drifting_gratings_75_repeats__contrast0-1';%'Signal_Averaged_dot_motion__Speed0-01--------0-02--------0-04'; %
FileName = 'Drifting_Grating';

load(fullfile(Path,File1),'Signal_FFT_Averaged');
Signal_FFT_Averaged_C08 = Signal_FFT_Averaged;

load(fullfile(Path,File2),'Signal_FFT_Averaged');
Signal_FFT_Averaged_C01 = Signal_FFT_Averaged;

clear Signal_FFT_Averaged;

%% load plotting data

load ROInames.mat;

ROIs = fieldnames(Signal_FFT_Averaged_C01);
ROIs = ROIs([4 3 6 1 5 2]);

COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs,'SubColors');
%% plot overlapped power spectrum phase locked

% you can plot this figure in two modes, non-normalized (MODE=1) or
% base-line removed (MODE=2)

MODE=1;

FIG = figure;
set(FIG,'color','w','unit','inch','position',[0 0 20 15])


for roi = 1:numel(ROIs)
    S1      = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMEvoked),3));
    STD1 = std(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMEvoked),[],3)./sqrt(size(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase,3));
    S1Pres  = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase),3));
    STD1Pres  = std(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase),[],3)./sqrt(size(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase,3));
    
    S8      = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMEvoked),3));
    STD8 = std(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMEvoked),[],3)./sqrt(size(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase,3));
    S8Pres  = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMPresPhase),3));
    STD8Pres  = std(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMPresPhase),[],3)./sqrt(size(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase,3));
    
    F       = Signal_FFT_Averaged_C08.(ROIs{roi}).Freq;
    FPres   = Signal_FFT_Averaged_C08.(ROIs{roi}).FreqPres;
    for l = 1:6
        subplot(6,numel(ROIs),(l-1)*numel(ROIs)+roi);
        if MODE==2
            LG(1) = plot(F,(squeeze(S1(l,:))-squeeze(S1Pres(l,:))),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5);
            hold on; LG(2) = plot(F,(squeeze(S8(l,:))-squeeze(S8Pres(l,:))),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5);
            ylim([0 max(S8(:)-S8Pres(:))])
        elseif MODE==1
            hold on;
            fill([F flip(F)],[S1(l,:)+STD1(l,:) flip(S1(l,:))-flip(STD1(l,:))],squeeze(Colors(roi,4,:))','facealpha',.3,'edgecolor','none')
            fill([F flip(F)],[S8(l,:)+STD8(l,:) flip(S8(l,:))-flip(STD8(l,:))],squeeze(Colors(roi,4,:))'./(2.5),'facealpha',.3,'edgecolor','none')
            fill([F flip(F)],[S1Pres(l,:)+STD1Pres(l,:) flip(S1Pres(l,:))-flip(STD1Pres(l,:))],'k','facealpha',.3,'edgecolor','none')
            LG(1)=plot(F,(S1(l,:)),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5);
            LG(2)=plot(F,squeeze(S8(l,:)),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5);
            LG(3)=plot(FPres,squeeze(S1Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5);
            %plot(FPres,squeeze(S8Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
            ylim([0 max(S8(:))])
        end
        xlim([1 50])
        
        if l==1
            title(ROI_names.(ROIs{roi}));
            legend(LG,{'Contrast=.1','Contrast=.8' ,'Prestimulus'})
            %legend({'Fast','Slow' ,'Prestimulus'})
        end
        if roi==1
            ylabel(['L' num2str(l)]);
        end
        if (roi==numel(ROIs)) && (l==6)
            
            xlabel('Frequency (Hz)');
        end
    end
    
end
if MODE==1
	export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Evoked_' FileName '_comparison']),'-pdf','-r200');
elseif MODE==2
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Evoked_' FileName '_comparison_BSremoved']),'-pdf','-r200');
end
%% plot overlapped power spectrum spontaneous

% you can plot this figure in two modes, non-normalized (MODE=1) or
% base-line removed (MODE=2)

MODE=2;
FIG = figure;
set(FIG,'color','w','unit','inch','position',[0 0 20 15])

for roi = 1:numel(ROIs)
    clear LG;
    Data1       = abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMSpon);
    Data1Pres   = abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPres);
    if MODE==2
        Data1 = Data1-Data1Pres;
    end
        
    S1      = abs(mean(Data1,3));
    S1Pres  = abs(mean(Data1Pres,3));
    STD1    = std(Data1,[],3)./sqrt(size(Data1,3));
    STD1Pres  = std(Data1Pres,[],3)./sqrt(size(Data1Pres,3));

    Data8       = abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMSpon);
    Data8Pres   = abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMPres);
    if MODE==2
        Data8   = Data8-Data8Pres;
    end
    S8      = abs(mean(Data8,3));
    S8Pres  = abs(mean(Data8,3));
    STD8    = std(Data8,[],3)./sqrt(size(Data8,3));
    STD8Pres  = std(Data8Pres,[],3)./sqrt(size(Data8Pres,3));

    F       = Signal_FFT_Averaged_C08.(ROIs{roi}).Freq;
    FPres   = Signal_FFT_Averaged_C08.(ROIs{roi}).FreqPres;
    for l = 1:6
        subplot(6,numel(ROIs),(l-1)*numel(ROIs)+roi);
        if MODE==2
            hold on;           
            fill([F flip(F)],[S1(l,:)+STD1(l,:) flip(S1(l,:))-flip(STD1(l,:))],squeeze(Colors(roi,4,:))','facealpha',.3,'edgecolor','none')
            fill([F flip(F)],[S8(l,:)+STD8(l,:) flip(S8(l,:))-flip(STD8(l,:))],squeeze(Colors(roi,4,:))'./(2.5),'facealpha',.3,'edgecolor','none')
            LG(1) = plot(F,squeeze(S1(l,:)),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5);
            LG(2) = plot(F,squeeze(S8(l,:)),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5);
            ylim([min(S8(:)) max(S8(:))])
            hline(0,'k--')
        elseif MODE==1
            hold on;
            fill([F flip(F)],[S1(l,:)+STD1(l,:) flip(S1(l,:))-flip(STD1(l,:))],squeeze(Colors(roi,4,:))','facealpha',.3,'edgecolor','none')
            fill([F flip(F)],[S8(l,:)+STD8(l,:) flip(S8(l,:))-flip(STD8(l,:))],squeeze(Colors(roi,4,:))'./(2.5),'facealpha',.3,'edgecolor','none')
            fill([F flip(F)],[S1Pres(l,:)+STD1Pres(l,:) flip(S1Pres(l,:))-flip(STD1Pres(l,:))],'k','facealpha',.3,'edgecolor','none')
            LG(1)=plot(F,(S1(l,:)),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5);
            LG(2)=plot(F,squeeze(S8(l,:)),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5);
            LG(3)=plot(FPres,squeeze(S1Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5);
            %plot(FPres,squeeze(S8Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
            ylim([0 max(S8(:))])
        end
        xlim([1 100])
        ylim([0 max(S8Pres(:))]*1.2)
        if l==1
            title(ROI_names.(ROIs{roi}));
            legend(LG,{'Contrast=.1','Contrast=.8' ,'Prestimulus'})
            %legend({'Fast','Slow' ,'Prestimulus'})
        end
        if roi==1
            ylabel(['L' num2str(l)]);
        end
        if (roi==numel(ROIs)) && (l==6)
            
            xlabel('Frequency (Hz)');
        end
    end
    
end
if MODE==1
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Induced_' FileName '_comparison']),'-pdf','-r200');
elseif MODE==2
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs_ASD_Induced_' FileName '_comparison_BSremoved']),'-pdf','-r200');
end

%%
