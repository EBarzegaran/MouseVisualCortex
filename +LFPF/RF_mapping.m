function [RF_mapping] = RF_mapping(gfilepath)

if exist(gfilepath,'file')
    Gabors = load(gfilepath);
else
    error('First import gabor data from AllenSDK');
end

%%
Pnames = {'x_position','y_position','orientation'};


XVals = sort(unique(Gabors.cnd_info.x_position));
YVals = sort(unique(Gabors.cnd_info.y_position));
OVals = sort(unique(Gabors.cnd_info.orientation));
TimeWin = Gabors.Times>.0 & Gabors.Times<.12;
for x = 1:numel(XVals)
    StimParams.x_position = XVals(x);
    for y = 1:numel(YVals)
        StimParams.y_position = YVals(y);
        for O = 1:numel(OVals)
            StimParams.orientation = OVals(O);
            CondSel = cellfun(@(x) ismember(Gabors.cnd_info.(x),StimParams.(x)),Pnames,'uni',false);
            CondSel = prod(cat(3,CondSel{:}),3);
            RF(x,y,O,:,:) = sqrt(squeeze(mean(Gabors.Y(:,:,TimeWin,CondSel==1).^2,3)));
        end
    end
end

%%  (1) first indicate the RF significance -> permutation
% RF is a X x Y x Orientation x trial x Electrodes
for dir = 1:3
    for E = 1:size(RF,5)
        ORF = squeeze(mean(RF(:,:,dir,:,E),4));
        OrStat = RFstat(ORF);
        RandDis =[];
        for p = 1:1000
            for t = 1:size(RF,4)
                ind = randperm(numel(ORF));
                RF_temp = RF(:,:,dir,t,E);
                RF_rand(:,:,t) = reshape(RF_temp(ind),size(ORF,1),size(ORF,2));
            end
            RandDis(p) = RFstat(mean(RF_rand,3));
        end
        Pval(dir,E) = sum(RandDis>OrStat)./numel(RandDis);
    end
end

%% (2) for significant electrodes calculate Rf size, and center and map
RFM = squeeze(mean(RF,4));

for dir = 1:3
    % FIG = figure;
    %set(FIG,'unit','inch','position',[1 1 5 15])
    for E = 1:size(RFM,4)
        IM = imgaussfilt(squeeze(RFM(:,:,dir,E)),.8);
        IM = IM./mean(IM);
        
%         if E<80
%             subplot(16,5,E),
%             SB = imagesc(IM);
%             Transp = 1-mean(Pval(:,E));
%             Transp(Transp<.8)=0;
%             set(SB,'alphadata',Transp)
%         end
        
        % Threshold
        TH(E,dir) = max(IM(:))-(std(IM(:))*1.2);
        
        %Area 
        CC = bwconncomp(IM>TH(E,dir));
        CS = cellfun(@(x) numel(x), CC.PixelIdxList);
        [Area(dir,E), Ind] = max(CS);
        
        IM2 = zeros(9,9);
        IM2(CC.PixelIdxList{Ind}) = IM(CC.PixelIdxList{Ind});
        
        % centroids
        props = regionprops(IM2>0, IM2, 'WeightedCentroid');
        RF_centroid(:,dir,E) = props.WeightedCentroid;
          
%         if E<80
%             %subplot(16,5,E),imagesc(IM2);
%         end
    end

end
%% return values
% CHECK THE DIMENSIONS
RF_mapping.Area     = Area;
RF_mapping.Centroid = RF_centroid;
RF_mapping.PValue   = Pval;
RF_mapping.Maps     = RFM;
RF_mapping.Conds.Xvals    = XVals;
RF_mapping.Conds.Yvals    = YVals;
RF_mapping.Conds.Ovals    = OVals;
end

function Stat = RFstat(MAP)
    E = mean(MAP);
    Stat = sum((MAP-E).^2)/E.^2;
end

