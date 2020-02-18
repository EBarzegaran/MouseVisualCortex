function PDC_boot = bootstrap_PDC (PDC,nboots)
% bootstraps a within area PDC matrix, to extract the significant
% connections

% PDC: is a node x node x frequency x time x subjects matrix
%% default variables

if ~exist('nboots','var') || isempty(nboots)
    nboots = 500;
end

%% prepare bootstrapping

Dims = size(PDC);
boots = randi(Dims(5),Dims(5),nboots);

for p = 1:nboots
    PDC_boot(:,:,:,:,p) = mean(PDC(:,:,:,:,boots(:,p)),5);
end

end