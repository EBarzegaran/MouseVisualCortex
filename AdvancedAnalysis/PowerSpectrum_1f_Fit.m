%% check the power spectrum maybe?
%Question why nlinfit does not work?
FileNames2 = {'Signal_PSD_drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15.mat',...
    'Signal_PSD_drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15.mat'};
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%
close all
model_fun.power_law = @(p,x) ( (p(2)*x+1).^-p(1)) ;
%model_fun.power_law = @(p,x) ( (x).^-p(1)) ;

TW = [-.3 0; .05 .5];
Colors = {'r','b'}
for i = 1:2
    PSD{i} = load(fullfile(Path,FileNames2{i}));
    figure
    for T = 1:2
        
        t = PSD{i}.t;
        Freqs = PSD{i}.Freqs;
        P_temp = PSD{i}.STOKPSD.VISp;
        P_temp = cat(4,P_temp{:});
        PP_temp = mean(P_temp(:,:,t<TW(T,2) & t>TW(T,1)),3);
        P = mean(PP_temp);

        %P = P-mean(P);
        opts = statset('nlinfit');
        opts.RobustWgtFun = 'bisquare';
        [this_p,R,J,CovB,MSE] =nlinfit(Freqs(4:75),P(4:75),model_fun.power_law,[3.,1],opts);
        param(T) = this_p(1)

        scatter(log10(Freqs(4:100)),log10(mean(PP_temp(:,4:100))),[],Colors{T})
        hold on;
        y = (this_p(2)*Freqs(4:100)+1).^-(this_p(1)) ;
        %y = (Freqs(4:100)).^-this_p(1) ;
        PL(T) = plot(log10(Freqs(4:100)),log10(y),Colors{T});
    end
    legend(PL,['Pre = ' num2str(round(param(1),2))],['Post = ' num2str(round(param(2),2))])
    ylim([-21 -11])
end
%%
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098';%'_dot_motion__Speed0-01--------0-02--------0-04_iPDC_Mord15';%'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord10';%
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/dynet_toolbox-master'));
fileparts(mfilename('fullpath'))

load([Path 'STOK_Average_' FileName '.mat']);
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
SavePath = Path;
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs);
load ROInames;
NROIs = numel(STOK_avg.ROIs);
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
% plot each ROI separately


Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;


%% Only for prestimulus time window
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
FigPath = fullfile(Path,'StatResults','PARAFAC');

clear this_p MSE
close all
FIG = figure;
set(FIG,'unit','inch','position',[5 5 8 6],'color','w')
FS = 12;
model_fun.power_law = @(p,x) ( p(2)*(x).^-p(1)) ;
TW = [-.3 0; .05 .5];
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
for roi1 = 1:NROIs
    ind1 = (roi1-1)*6+1:roi1*6;    
    for roi2 = 1:NROIs
        clear pdct
        if roi1~= roi2 % if PDC
            ind2 = (roi2-1)*6+1:roi2*6;
            pdct = reshape(PDC(ind1,ind2,:,:),36,size(PDC,3),size(PDC,4));
            pdct = mean(pdct(:,:,Time>TW(1,1) & Time<TW(1,2)),3);
            pdct = (pdct-min(pdct,[],2))./(max(pdct,[],2)-min(pdct,[],2));
            %pdct2 = interp1(Freq,pdct',.5:.5:100);
            %[this_p,~,~,~,MSE] =arrayfun(@(x) nlinfit(Freq,mean(pdct),model_fun.power_law,[3.,1],opts),1:36,'uni',false);
            [this_p,~,~,~,MSE] =nlinfit(Freq,nanmean(pdct),model_fun.power_law,[3.,1],opts);
            %param   = cat(1,this_p{:});
            %MSE     = cat(1,MSE{:});
            
            subtightplot(NROIs,NROIs,roi1+(roi2-1)*NROIs);
            plot(log10(Freq),nanmean(log10(pdct)),'color',[.2 .2 .2],'linewidth',1.5)
            
            y = this_p(2)*(Freq).^-(this_p(1)) ;
            hold on;
            plot(log10(Freq),log10(y),'--r','linewidth',1.5)
            text(0.2,-1.7,['\beta=' num2str(round(this_p(1),1))],'color','k')
            %axis tight;
            ylim([-2.2 .5])
            if roi2~=NROIs || roi1~=1
                set(gca,'xticklabels',[],'yticklabels',[]);
            else
                Fs = [1 10 100];
                set(gca,'xtick',log10(Fs),'xticklabels',Fs);
            end
        else % if PSD, diagonal
            for i = 1:numel(ind1)
                pdct(i,:,:) = PDC(ind1(i),ind1(i),:,:);
            end
            pdct = mean(pdct(:,:,Time>TW(1,1) & Time<TW(1,2)),3);
            pdct = (pdct-min(pdct,[],2))./(max(pdct,[],2)-min(pdct,[],2));
            [this_p,~,~,~,MSE] =nlinfit(Freq,nanmean(pdct),model_fun.power_law,[3.,1],opts);


            subtightplot(NROIs,NROIs,roi1+(roi2-1)*NROIs);
            plot(log10(Freq),nanmean(log10(pdct)),'color','b','linewidth',1.5)
            y = this_p(2)*(Freq).^-(this_p(1)) ;
            hold on;
            plot(log10(Freq),log10(y),'--r','linewidth',1.5)
            text(0.2,-4.5,['\beta=' num2str(round(this_p(1),1))],'color','b')
            %axis tight;
            ylim([-5 .5])
            if roi1~=1
                set(gca,'xticklabels',[],'yticklabels',[])
            else
                
                 set(gca,'xticklabels',[]);
            end
            if roi1==NROIs
                text(2.2,-5.4,'log10(PSD)   log10(iPDC)','rotation',90,'fontsize',FS)
                xlabel('Frequency(Hz)');
            end
            
        end
        set(gca,'linewidth',1.5,'fontsize',FS)
        if roi2==1
            title(ROISN{roi1})
        end
        if roi1==1
            YL = ylabel(ROISN{roi2},'fontweight','bold');
            if roi2==1
                PYL = get(YL,'position');
            end
            YL.Position(1)=PYL(1)+.15;
        end
    end
end


export_fig(FIG,fullfile(FigPath,[FileName '_Beta_Exponent']),'-pdf','-r200')

