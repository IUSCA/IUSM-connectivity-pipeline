function f_get_largest_clusters_only(fileIn,fileOut,threshold)

v=MRIread(fileIn);

N = max(v.vol(:));
vol_clean = zeros(size(v.vol));

for i=1:N
    clusters = bwconncomp(v.vol==i);
    for j=1:clusters.NumObjects
        if length(clusters.PixelIdxList{j})>threshold
            vol_clean(clusters.PixelIdxList{j}) = i;
        end
    end
end

v.vol = vol_clean;
MRIwrite(v,fileOut);

