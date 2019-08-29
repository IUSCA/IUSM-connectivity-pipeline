% ts are the time series
% numComps is the number of PCA components desired
function [PCtop,variance] = get_pca(ts,numComps)


if numComps>0
    [~,PC,latent] = pca(ts);
    variance = cumsum(latent)./sum(latent); %explained variance
    PCtop = PC(:,1:numComps);
    variance = variance(1:numComps);
else
    PCtop = [];
    variance = [];
end