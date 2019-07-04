% function resid = apply_regressors(ts,regressors)
% 
% N = size(ts,1);
% resid = zeros(size(data));
% 
% for i=1:N
%     currentTS = ts(i,:);
%     [~,~,currentResid] = regress(currentTS,regressors);
%     resid(i,:) = currentResid;
% end

function resid = apply_regressors(data,mask,regressors,scrubbing)


if nargin<4
    scrubbing = true(size(data,4),1);
end

regressors = unique(regressors','rows')'; %remove identical regressors if present

[sizeX,sizeY,sizeZ,numTimePoints] = size(data);
resid = zeros(sizeX,sizeY,sizeZ,numTimePoints);

for i=1:sizeX
    for j=1:sizeY
        for k=1:sizeZ
            if mask(i,j,k)
                TSvoxel = reshape(data(i,j,k,:),[numTimePoints,1]);
                B = regress(TSvoxel(scrubbing),regressors(scrubbing,:)); % coeffs learned from good points only
                B = repmat(B,1,numTimePoints);
                Yhat = sum(B.*regressors');
                resid(i,j,k,:) = TSvoxel-Yhat';
            end
        end
    end
    if (mod(i,25)==0)
        disp(i/sizeX)
    end
end