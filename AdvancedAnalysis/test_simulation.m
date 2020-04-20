clear; clc;
close all
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'))
rmpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim/External/eeglab14_1_1b'))
FigPath = '';

%% design network
Net1 = BrainNetSim(3,100,.96);
Net1 = Net1.AddNode('AR');
% define net dynamics
NS = 100; % length of simulation: time points

Net1 = Net1.AddNodeFreqs([1 2 3],{[9],[4 15],[5 12]});
Net1 = Net1.GenerateARMatrix;
% [Net1,TS_un] = Net1.Realization(NS);
% TS_un = TS_un./sqrt(sum(TS_un.^2,2))*50;
for tr =1:100
    [~,TS_un(tr,:,:)] = Net1.Realization(NS);
end
TS_un = TS_un./sqrt(sum(TS_un.^2,3))*20;

plot_net_sig(Net1,reshape(permute(TS_un,[2 3 1]),size(TS_un,2),size(TS_un,1)*size(TS_un,3)),1,FigPath,'Simulation_temp_un');
% estimate the PDC
KF = dynet_SSM_KF(TS_un(:,:,:),6,.9);

S_STOK      = dynet_ar2pdc(KF,100,1:20,'iPDCs',1, 2,0);

dynet_connplot(S_STOK,1:200,1:20,[],[],[],[],1)

export_fig('Simulation_STOK_un','-pdf')

%%
% Add connections to the network
Net1 = Net1.AddConnection([1 2],'Type','bandpass','LF',5,'HF',15,'Order',6,'Gain',1/5);
Net1 = Net1.AddConnection([2 3],'Type','bandpass','LF',5,'Order',6,'Gain',1/5);
Net1 = Net1.AddConnection([3 4],'Type','high','HF',7,'Gain',1/5);
%Net1 = Net1.AddConnection([1 4],'Type','delay','Order',5,'Gain',1/5);

Net1 = Net1.GenerateARMatrix;
% Net1.ARMatrix(2,2,:) = Net1.ARMatrix(2,3,:);
% Net1.Nodes(2).Alpha = 0.;
% Net1.Nodes(3).Alpha = 0.;
% Net1.Nodes(4).Alpha = 0.;
clear TS
for tr =1:100
    tr
    [~,TS(tr,:,:)] = Net1.Realization(NS);
end

TS = TS./sqrt(sum(TS.^2,3))*20;
%
plot_net_sig(Net1,reshape(permute(TS,[2 3 1]),size(TS,2),size(TS,1)*size(TS,3)),1,FigPath,'Simulation_temp');

%% estimate the PDC
KF = dynet_SSM_KF(TS(:,1:4,:),6,.9);

S_STOK      = dynet_ar2pdc(KF,100,1:20,'iPDCs',1, 2,0);

KF_test = KF;
KF_test.AR = permute(repmat(Net1.ARMatrix,[1 1 1 10]),[2 1 3 4]);
S_STOK_test      = dynet_ar2pdc(KF_test,100,1:20,'iPDCs',1, 2,0);

dynet_connplot(S_STOK,1:200,1:20,[],[],[],[],1)
dynet_connplot(S_STOK_test,1:200,1:20,[],[],[],[],1)

export_fig('Simulation_STOK','-pdf')

%% extra plots
for i = 1:4
    subplot(4,1,i)
    PSD = sqrt(squeeze(S_STOK(i,i,:,20:end)));
    imagesc(PSD);%./sum(S_STOK(i,:,:,20:end),2)));
    caxis([0 max(PSD(:))/100]);
    axis xy
end

%% Maybe add 1/f noise
VTS = mean(std(TS,[],3));
WN = rand(size(TS)).*VTS;
WNFFT = fft(WN,[],3);
func = 1./(1:ceil((size(WNFFT,3)-1)/2));
PNFFT = WNFFT.* permute(repmat([0 func flip(func(1:end-1))]',[1 4 100]),[3 2 1]);
PN = ifft(PNFFT,[],3);


%%
figure
for i =1:4  
    PSD = squeeze(sum(S_STOK(:,i,:,:)));
    subplot(4,1,i),
    imagesc(PSD);
    axis xy
    %caxis([0 max(PSD(:))/10000])
end

function plot_net_sig(Net,TS,issave, FigPath,Name)
    % Plots Network realization
    if ~exist('issave','var') || isempty(issave)
        issave = 0;
    end
    if ~exist('FigPath','var') || isempty(FigPath)
        FigPath = pwd;
    end
    if ~exist('Name','var') || isempty(Name)
        Name = 'NetResults';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fig1 = figure;
    for node = 1:Net.NodeNum
        % 
        subplot(Net.NodeNum,2,(node-1)*2+1),plot(TS(node,1:Net.SF*2));
        set(gca,'xtick',0:Net.SF:Net.SF*2,'xticklabel',0:2)
        if node == Net.NodeNum, xlabel('time (s)'); end
        ylabel(['Node' num2str(node)],'fontweight','bold','fontsize',11);
        if node ==1, title('Temporal dynamic'); end

        subplot(Net.NodeNum,2,(node-1)*2+2);
        [Z,f] = pwelch(TS(node,:),Net.SF,[],[],Net.SF);
        plot(f,Z,'linewidth',2);
        if node==1, title('Power Spectrum Density');end
        if node == Net.NodeNum,  xlabel('Frequency(Hz)');end
        xlim([0 20]);
        ylim([0 1])
        set(gca,'yticklabel',[])
        %legend(['Signal #' num2str(Sig)])
    end
    if issave
        set(Fig1,'PaperPosition',[1 1 6 4.5])
        print(fullfile(FigPath,Name),'-r300','-dtiff')
    end
end

