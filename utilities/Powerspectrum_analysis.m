clc;

Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged';

load(fullfile(Path,'Fullmodel','Signal_PSD'));
%%
for S = 1:1
    figure(1);
    for roi = 1:numel(ROIs)
        PSD_temp = WaveletPSD.(ROIs{roi}){S};        
        if ~isempty(PSD_temp)
            PSD_temp = (PSD_temp - mean(PSD_temp(:,:,t<0 & t<-.3),3))./mean(PSD_temp(:,:,t<0 & t<-.3),3);
            for l =1:6
                subplot(numel(ROIs),6,(l-1)*6+roi)
                
                imagesc(squeeze(PSD_temp(l,:,:)));
                axis xy;
                caxis([0 4])
                if l==1, title(ROIs{roi}); end
                %caxis([0 10.^-8]);
                colormap('jet');
            end
        end
    end
end
