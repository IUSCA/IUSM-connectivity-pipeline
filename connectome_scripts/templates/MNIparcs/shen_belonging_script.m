clearvars
shen= MRIread('shen_MNI152_v1_5.nii.gz');
yeo7=MRIread('yeo7_MNI152.nii.gz');
yeo17=MRIread('yeo17_MNI152.nii.gz');

shenVOL=shen.vol;
yeo7VOL=yeo7.vol;
yeo17VOL=yeo17.vol;
% ----------------------------------------------------%
shenSizes=NaN(286,1);
yeo7Sizes=NaN(7,1);
yeo17Sizes=NaN(17,1);

shenYEO7vox=NaN(286,7);     % HARD CODED SIZES!!!!
shenYEO7per=NaN(286,7);
shenYEO7leftoversVOX=NaN(286,1);
shenYEO7leftoversPER=NaN(286,1);

shenYEO17vox=NaN(286,17);
shenYEO17per=NaN(286,17);
shenYEO17leftoversVOX=NaN(286,1);
shenYEO17leftoversPER=NaN(286,1);
% ----------------------------------------------------%
for i=1:max(max(max(shenVOL)))
    rMask=shenVOL==i;
    rSize=nnz(rMask);
    shenSizes(i,1)=rSize;
    
    for j=1:max(max(max(yeo7VOL)))
        if i==1   
            nMask=yeo7VOL==j;
            nSize=nnz(nMask);
            yeo7Sizes(j,1)=nSize;
        end
    end
    ShenYeos=yeo7VOL(rMask);
    Yeonique=unique(ShenYeos);
    for k=1:size(Yeonique,1)
        tvox=nnz(ShenYeos==Yeonique(k,1));
        tper=(tvox/rSize)*100;
        if Yeonique(k,1)==0
            shenYEO7leftoversVOX(i,1)=tvox;
            shenYEO7leftoversPER(i,1)=tper;
        else
            shenYEO7vox(i,Yeonique(k,1))=tvox;
            shenYEO7per(i,Yeonique(k,1))=tper;
        end
     end
    
    for j=1:max(max(max(yeo17VOL)))
        if i==1   
            nMask=yeo17VOL==j;
            nSize=nnz(nMask);
            yeo17Sizes(j,1)=nSize;
        end
    end
    ShenYeos=yeo17VOL(rMask);
    Yeonique=unique(ShenYeos);
    for k=1:size(Yeonique,1)
        tvox=nnz(ShenYeos==Yeonique(k,1));
        tper=(tvox/rSize)*100;
        if Yeonique(k,1)==0
            shenYEO17leftoversVOX(i,1)=tvox;
            shenYEO17leftoversPER(i,1)=tper;
        else
            shenYEO17vox(i,Yeonique(k,1))=tvox;
            shenYEO17per(i,Yeonique(k,1))=tper;
        end
    end
end

save('shen_1_5_belonging_to_yeo_MNI152.mat','shenSizes','yeo7Sizes','yeo17Sizes',...
    'shenYEO7vox','shenYEO7per','shenYEO7leftoversVOX','shenYEO7leftoversPER',...
    'shenYEO17vox','shenYEO17per','shenYEO17leftoversVOX','shenYEO17leftoversPER');

        