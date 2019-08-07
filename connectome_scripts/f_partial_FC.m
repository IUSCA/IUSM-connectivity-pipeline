function FCpartialFull = f_partial_FC(restingROIs,validROIs,scrubbing,noiseTerm)

if nargin<3
    scrubbing = false(size(restingROIs,2),1);
end
if nargin<4
    noiseTerm = 0.001;
end
Nfull = size(restingROIs,1);
N = nnz(validROIs);
FCaux = inv(cov(restingROIs(validROIs,~scrubbing)')+(eye(N,N)*noiseTerm));
FCpartial = zeros(size(FCaux));

for i=1:N
    for j=i+1:N
        if i~=j
            FCpartial(i,j) = -FCaux(i,j)/sqrt(FCaux(i,i)*FCaux(j,j));
        end
    end
end

FCpartial = FCpartial + FCpartial';

FCpartialFull = nan(Nfull,Nfull);
FCpartialFull(validROIs,validROIs)=FCpartial;


        