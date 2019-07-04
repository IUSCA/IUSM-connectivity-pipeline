%% Parcellation Region of Interest (ROI) belonging

%   For a given input parc1 (nifti volume with integer ROI labels), outputs
%   what portion of each ROI belongs to integer labels of parc2.
%   -Evgeny Chumin, Indiana University School of Medicine, 2019

%   REQUIRED: MRIread.m - Part of the FreeSurfer Package.
%------------------------------------------------------------------------%
 clearvars
%---------%
%% INPUTS
parc1= MRIread('shen_MNI152_v1_5.nii.gz');
parc2=MRIread('yeo7_MNI152.nii.gz');
%-----------------------------------------%
% determine parcellation sizes
vol1=parc1.vol;
vol2=parc2.vol;
mxparc1=max(max(max(vol1)));
mxparc2=max(max(max(vol2)));

%% OUTPUTS
% Size of all ROI in parc1 and 2
parc1Sizes=NaN(mxparc1,1);
parc2Sizes=NaN(mxparc2,1);

% number and percent of voxels of ROI in parc1 that overlap ROI of parc2
parc1_2_vox=NaN(mxparc1,mxparc2);
parc1_2_per=NaN(mxparc1,mxparc2);

% number of voxels of ROI in parc1 that did not match an index in parc2
parc1_2_leftoversVOX=NaN(mxparc1,1);
parc1_2_leftoversPER=NaN(mxparc1,1);

% ----------------------------------------------------%
for i=1:mxparc1 % for every parc1 region
    rMask=vol1==i; % make logic mask
    rSize=nnz(rMask); % get size of ROI in voxels
    parc1Sizes(i,1)=rSize; 
    
    % get size of all parc2 indices
    for j=1:mxparc2 
        if i==1   
            nMask=vol2==j;
            nSize=nnz(nMask);
            parc2Sizes(j,1)=nSize;
        end
    end
    %------------------------------%
    roiIdx=vol2(rMask); % get parc2 indices of parc1 region
    roi_unique=unique(roiIdx);
    for k=1:size(roi_unique,1) % for every parc2 label in parc1 roi
        tvox=nnz(roiIdx==roi_unique(k,1)); % # of voxels in roi with parc2 label k
        tper=(tvox/rSize)*100; % percent of voxels
        % if parc2 idx=0, bin into leftovers
        if roi_unique(k,1)==0
            parc1_2_leftoversVOX(i,1)=tvox;
            parc1_2_leftoversPER(i,1)=tper;
        else % else plase it into ([row]parc1,[col]parc2)
            parc1_2_vox(i,roi_unique(k,1))=tvox;
            parc1_2_per(i,roi_unique(k,1))=tper;
        end
    end
end
clear i j k mxparc1 mxparc2 parc1 parc2 vol1 vol2 nSize rSize nMask rMask ...
    tper tvox roiIdx roi_unique


        