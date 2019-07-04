function bvecs_corr = f_correct_bvec_from_matrix(bvecs,matrix)

if size(bvecs,1)==3
    bvecs = bvecs';
end

numDirs = size(bvecs,1);
bvecs_corr = nan(numDirs,3);
for dwi_dir=1:numDirs
    bvecs_corr(dwi_dir,:) = matrix(1:3,1:3) * bvecs(dwi_dir,:)';
end