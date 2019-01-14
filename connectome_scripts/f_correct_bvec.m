function [bvecs_corr,angular_motion_deg,angular_motion_rad] = f_correct_bvec(fileRef,fileDWI,file_bvec,fileDWI_corr,file_bvec_corr,path2DWI,FSLsetup,path2FSL)

    bvecs = dlmread(file_bvec,',');
    numVols = max(size(bvecs));
    
    if size(bvecs,1)==3
        bvecs = bvecs';
    end

    bvecs_corr = zeros(size(bvecs));
    
     path2split = fullfile(path2DWI,'3_DWI_Split');
    if ~exist(path2split,'dir')
        mkdir(path2split);
    end
    fileMatrix = fullfile(path2split,'temp_matrix.txt');
    filesCorr = ''; % initialization
    sentence = sprintf('%s;%s/fslsplit %s %s/',FSLsetup,path2FSL,fileDWI,path2split);
    [status,result] = system(sentence);
    
    for dwi_dir = 1:numVols
        dwi_dir
        
    %% obtain fileNames
        if dwi_dir<11
            suffix = sprintf('000%d',dwi_dir-1);
        elseif dwi_dir<101
            suffix = sprintf('00%d',dwi_dir-1);
        else
            suffix = sprintf('0%d',dwi_dir-1);
        end
        
        fileSingleVol = fullfile(path2split,sprintf('%s.nii.gz',suffix)); % original volume of a single diffusion direction
        fileSingleVolCorr = fullfile(path2split,sprintf('%s_corr.nii.gz',suffix)); % volume registered to b0 (rigid body)
        
        %% rigid body registration to b0
        sentence = sprintf('%s;%s/flirt -in %s -ref %s -nosearch -o %s -dof 6 -omat %s -interp spline',...
        FSLsetup,path2FSL,fileSingleVol,fileRef,fileSingleVolCorr,fileMatrix);
        [status,result] = system(sentence);
        
        %% update bvec
        matrix = dlmread(fileMatrix);
        matrix(1:3,1:3)
        bvecs_corr(dwi_dir,:) = matrix(1:3,1:3) * bvecs(dwi_dir,:)';
        
        %% delete registration matrix
        clear matrix;
        sentence = sprintf('rm %s',fileMatrix);
        [status,result] = system(sentence);
        %% append vol to be added with fslmerge
        filesCorr = sprintf('%s %s',filesCorr,fileSingleVolCorr);
            
    end
    
    sentence = sprintf('%s;%s/fslmerge -t %s %s ',FSLsetup,path2FSL,fileDWI_corr,filesCorr);
    [status,result] = system(sentence);
    
    sentence = sprintf('rm -rf %s',path2split);
    [status,result] = system(sentence);
    
    dlmwrite(file_bvec_corr,bvecs_corr,'delimiter',',','precision','%.8f')
    
    angular_motion_rad = zeros(numVols,1);
    angular_motion_deg = zeros(numVols,1);
    for i=1:numVols
        dir_uncorr = bvecs(i,:);
        dir_corr = bvecs_corr(i,:);
        angular_motion_rad(i) = acos(sum(dir_uncorr.*dir_corr)./sqrt(sum(dir_uncorr.^2))./sqrt(sum(dir_corr.^2)));
        angular_motion_deg(i) = 180*angular_motion_rad(i)./pi;
    end