clear; 
%close all; 
clc;
Sessions_ID = '771990200';
P_ind       = 1;

ProjectPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/';
Stimuli     = '_dot_motion';
Probes_Stim = subfiles(fullfile(ProjectPath,Sessions_ID,'MatlabData',['*' Stimuli '250.mat']),0);
Probes_ID   = cellfun(@(x) x(1:subsref(strfind(x,'_'),struct('type','()','subs',{{1}}))-1),Probes_Stim,'uni',false);
Probes_ID   = unique(Probes_ID);


Motion  = load(['/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/' Sessions_ID '/MatlabData/' Probes_ID{P_ind} '_dot_motion250.mat']);
P_info  = load(['/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/' Sessions_ID '/MatlabData/' Probes_ID{P_ind} '_ProbeInfo.mat']);
Grating = load(['/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/' Sessions_ID '/MatlabData/' Probes_ID{P_ind} '_drifting_gratings_75_repeats250.mat']);

disp(cat(2,P_info.structure_acronyms{end-3:end-1}))
%% Dot motion
[Speed,Ind1] = sort(Motion.cnd_info.Speed);
[~,Ind2] = sort(Motion.cnd_info.Dir(Ind1));


Dir = Motion.cnd_info.Dir(Ind1);
Dir = Dir(Ind2);
 
 %Speeds = Motion.cnd_info.Speed(Ind1);
 Speed = Speed(Ind2);
 
 Y = Motion.Y(:,:,:,Ind1(Ind2));
 YN = Y - mean(Y(:,:,Motion.Times<0 & Motion.Times>-.25,:),3);
%  FIG= figure;
%  set(FIG,'unit','inch','position',[1 1 25 10])

%  for i = 1:28
%      subplot(4,7,i)
%      imagesc(squeeze(mean(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i)))./squeeze(std(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i))))
%      axis xy
%      %caxis([-1 1]*9*10^-5)
%      caxis([-1 1]*1.6)
%      if i <8
%          title(['Speed = ' num2str(Speed(i))]);
%      end
%      if mod(i,7)==1, ylabel([' - Dir = ' num2str(Dir(i))]);end
%      
%  end
% colormap(jmaColors('coolhot')) 
%%  speed tunning
FIG= figure;
set(FIG,'unit','inch','position',[1 1 10 15])

 Y = Motion.Y(:,:,:,Ind1(Ind2));
 YN = Y - mean(Y(:,:,Motion.Times<0 & Motion.Times>-.25,:),3);
 for i = 1:7
     subplot(7,1,i)
     imagesc(squeeze(mean(mean(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i:7:end)),4)));%./squeeze(std(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i))))
     %imagesc(squeeze(mean(mean(Y(:,:,:,i:7:end)),4)));
     axis xy
     caxis([-1 1]*2*10^-5)
     %caxis([-1 1]*1.6)
     if i <8
         title(['Speed = ' num2str(Speed(i))]);
     end
     
 end
colormap(jmaColors('coolhot')) 


%%  Orient tunning
if false
FIG= figure;
set(FIG,'unit','inch','position',[1 1 10 7])

 Y = Motion.Y(:,:,:,Ind1(Ind2));
 YN = Y - mean(Y(:,:,Motion.Times<0 & Motion.Times>-.25,:),3);
 for i = 1:4
     subplot(4,1,i)
     imagesc(squeeze(mean(mean(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,(1:7)+(7*(i-1)))),4)));%./squeeze(std(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i))))
     axis xy
     caxis([-1 1]*2*10^-5)
     %caxis([-1 1]*1.6)
     if i <8
         title(['Dir = ' num2str(Dir((i-1)*7+1))]);
     end
     vline(find(round(Motion.Times,2)==0,1),'w--')
 end
colormap(jmaColors('coolhot')) 
end
 %% Grating
if false
    FIG = figure;
    set(FIG,'unit','inch','position',[1 1 15 5])

    [Orient,Ind1] = sort(Grating.cnd_info.orientation);
    [~,Ind2] = sort(Grating.cnd_info.contrast(Ind1));


    Cont = Grating.cnd_info.contrast(Ind1);
    Cont = Cont(Ind2);

     %Speeds = Motion.cnd_info.Speed(Ind1);
    Orient = Orient(Ind2);

     Y = Grating.Y(:,:,:,Ind1(Ind2));
     YN = Y - mean(Y(:,:,Grating.Times<0 & Grating.Times>-.25,:),3); 
    for i = 1:8
         subplot(2,4,i)
         imagesc(squeeze(mean(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i)))./squeeze(std(YN(:,P_info.intervals(end-2):P_info.intervals(end-1),:,i))))
         axis xy
         %caxis([-1 1]*9*10^-5)
         caxis([-1 1])
         title(['Orient = ' num2str(Orient(i)) ' - Cont = ' num2str(Cont(i))])
         xlim([find(round(Grating.Times,2)==-.5,1) find(round(Grating.Times,2)==1,1)])
    end
    colormap(jmaColors('coolhot')) 
end
 
