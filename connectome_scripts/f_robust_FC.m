function [FCrobust,VProbust] = f_robust_FC(restingROIs,validROIs,scrubbing)

N = length(validROIs);

if size(restingROIs,2)==N
    restingROIs = restingROIs';
end

numPoints = size(restingROIs,2);

if nargin<3
    scrubbing = true(numPoints,1);
else
    scrubbing = logical(scrubbing);
end

FCrobust = zeros(N,N);
VProbust = zeros(N,N);


ts_outliers = f_detect_univariate_outliers(restingROIs,scrubbing');

ts_outliers = ~logical(ts_outliers');
restingROIs = restingROIs';

for i=1:N-1
    
    if validROIs(i)
        roi1 = restingROIs(:,i);
        validPoints_i = ts_outliers(:,i) & (scrubbing);
        for j=i+1:N
            if validROIs(i)
                roi2 = restingROIs(:,j);
                validPoints = validPoints_i & ts_outliers(:,j);
                if  nnz(validPoints)>1
                    FCrobust(i,j) = corr(roi1(validPoints),roi2(validPoints));
                end
                VProbust(i,j) = nnz(validPoints)/numPoints;
            end
        end
    end
end

FCrobust = FCrobust + FCrobust';
FCrobust(logical(eye(N,N))) = 0;
VProbust = VProbust + VProbust';
VProbust(logical(eye(N,N))) = 0;

        