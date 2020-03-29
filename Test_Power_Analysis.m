
clear; clc;

Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
File1 = 'Signal_Averaged_dot_motion__Speed0-0005-------0-001-------0-005.mat';% 'Signal_Averaged_drifting_gratings_75_repeats__contrast0-8'
File2 = 'Signal_Averaged_dot_motion__Speed0-01--------0-02--------0-04'; %'Signal_Averaged_drifting_gratings_75_repeats__contrast0-1'
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
%% plot overlapped power spectrum
FIG = figure;
set(FIG,'color','w','unit','inch','position',[0 0 20 15])

for roi = 1:numel(ROIs)
    S1      = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMEvoked),3));
    S1Pres  = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPresPhase),3));
    
    S8      = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMEvoked),3));
    S8Pres  = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMPresPhase),3));
    F       = Signal_FFT_Averaged_C08.(ROIs{roi}).Freq;
    FPres   = Signal_FFT_Averaged_C08.(ROIs{roi}).FreqPres;
    for l = 1:6
        subplot(6,numel(ROIs),(l-1)*numel(ROIs)+roi);
        plot(F,squeeze(S1(l,:)),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5)
        hold on; plot(F,squeeze(S8(l,:)),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5)
        
        plot(FPres,squeeze(S1Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
        %plot(FPres,squeeze(S8Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
        xlim([1 50])
        ylim([0 max(S8(:))])
        if l==1
            title(ROI_names.(ROIs{roi}));
            %legend({'Contrast=.1','Contrast=.8' ,'Prestimulus'})
            legend({'Fast','Slow' ,'Prestimulus'})
        end
        if roi==1
            ylabel(['L' num2str(l)]);
        end
        if (roi==numel(ROIs)) && (l==6)
            
            xlabel('Frequency (Hz)');
        end
    end
    
end

export_fig(FIG,fullfile(Path,'Signal_Averaged_AllROIs_ASD_Evoked_dot_motion_comparison'),'-pdf','-r200');

%% plot overlapped power spectrum
FIG = figure;
set(FIG,'color','w','unit','inch','position',[0 0 20 15])

for roi = 1:numel(ROIs)
    S1      = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMSpon),3));
    S1Pres  = abs(mean(abs(Signal_FFT_Averaged_C01.(ROIs{roi}).SMPres),3));
    
    S8      = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMSpon),3));
    S8Pres  = abs(mean(abs(Signal_FFT_Averaged_C08.(ROIs{roi}).SMPres),3));
    F       = Signal_FFT_Averaged_C08.(ROIs{roi}).Freq;
    FPres   = Signal_FFT_Averaged_C08.(ROIs{roi}).FreqPres;
    for l = 1:6
        subplot(6,numel(ROIs),(l-1)*numel(ROIs)+roi);
        plot(F,squeeze(S1(l,:)),'color',squeeze(Colors(roi,4,:)),'linewidth',1.5)
        hold on; plot(F,squeeze(S8(l,:)),'color',squeeze(Colors(roi,4,:))/2.5,'linewidth',1.5)
        
        plot(FPres,squeeze(S1Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
        %plot(FPres,squeeze(S8Pres(l,:)),'-.','color',[.3 .3 .3],'linewidth',1.5)
        xlim([1 50])
        ylim([0 max(S8Pres(:))])
        if l==1
            title(ROI_names.(ROIs{roi}));
            %legend({'Contrast=.1','Contrast=.8' ,'Prestimulus'})
            legend({'Fast','Slow' ,'Prestimulus'})
        end
        if roi==1
            ylabel(['L' num2str(l)]);
        end
        if (roi==numel(ROIs)) && (l==6)
            
            xlabel('Frequency (Hz)');
        end
    end
    
end

export_fig(FIG,fullfile(Path,'Signal_Averaged_AllROIs_ASD_Induced_dot_motion_comparison'),'-pdf','-r200');


