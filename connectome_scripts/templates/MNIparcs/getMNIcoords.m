C=[];
Cmni=[];
for roi=1:286
    aux = find(v.vol==roi);
    [i,j,k] = ind2sub(size(v.vol),aux);
    C(roi,:)=[mean(i),mean(j),mean(k)];
    aux = [round(C(roi,:)),1];
    Cmni(roi,:) = v.vox2ras*aux';
end
Cmni(:,4)=[];
