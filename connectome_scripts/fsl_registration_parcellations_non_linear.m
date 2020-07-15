function [paths,configs,parcs]=fsl_registration_parcellations_non_linear(paths,configs,parcs)
%              FSL_REGISTRATION_PARCELLATIONS_NON_LINEAR
% Transform subject T1 into standard space.
% Generate MNI <-> Native transformation matrices
% Transform parcellation volumes into subject native space.
%
% Contributors:
%   Joaquin Goni, Purdue University
%   Joey Contreras, University Southern California
%   Mario Dzemidzic, Indiana University School of Medicine
%   Evgeny Chumin, Indiana University School of Medicine

%%
% Set registration directory path and remove any existing reg directories.
paths.T1.reg = fullfile(paths.T1.dir,'registration');

% Check if existng transformation matrices should be used.
if configs.T1.useExistingMats == 1
    % check if transformations are there
    if exist(paths.T1.reg,'dir')
        fileMat_dof6 = fullfile(paths.T1.reg,'T12MNI_dof6.mat');
        if exist(fileMat_dof6,'file')
            fileMat_dof6_inv = fullfile(paths.T1.reg,'MNI2T1_dof6.mat');
            if exist(fileMat_dof6_inv,'file')
                fileMat_dof12 = fullfile(paths.T1.reg,'T12MNI_dof12.mat');
                if exist(fileMat_dof12,'file')
                    fileMat_dof12_inv = fullfile(paths.T1.reg,'MNI2T1_dof12.mat');
                    if exist(fileMat_dof12_inv,'file')
                        fileWarp = fullfile(paths.T1.reg,'T12MNI_warp.nii.gz');
                        if exist(fileWarp,'file')
                            fileWarpInv = fullfile(paths.T1.reg,'MNI2T1_warp.nii.gz');
                            if exist(fileWarpInv,'file')
                                disp('Using existing transformation matrices')
                            else
                                warning('%s missing: Running reg2MNI',fileWarpInv)
                                configs.T1.useExistingMats = 0;
                            end
                        else
                            warning('%s missing: Running reg2MNI',fileWarp)
                            configs.T1.useExistingMats = 0;
                        end
                    else
                        warning('%s missing: Running reg2MNI',fileMat_dof12_inv)
                        configs.T1.useExistingMats = 0;
                    end
                else
                    warning('%s missing: Running reg2MNI',fileMat_dof12)
                    configs.T1.useExistingMats = 0;
                end
            else
                warning('%s missing: Running reg2MNI',fileMat_dof6_inv)
                configs.T1.useExistingMats = 0;
            end
        else
            warning('%s missing: Running reg2MNI',fileMat_dof6)
            configs.T1.useExistingMats = 0;
        end
    else
        warning('%s missing: Running reg2MNI',paths.T1.reg)
        configs.T1.useExistingMats = 0;
    end
end % EJC 2018.03.21 Separated file checks above into their own loop
if configs.T1.useExistingMats == 0 || isempty(configs.T1.useExistingMats) % perform transformations
    
    if exist(paths.T1.reg,'dir')
        sentence = sprintf('rm -rf %s',paths.T1.reg);
        [~,result] = system(sentence);
        disp(result)
    end
    mkdir(paths.T1.reg)

%% Register T1 to MNI and obtain inverse transformations
%% flirt dof 6
    disp('     T1 --> MNI152')
    if configs.T1.useMNIbrain == 1
        fileRef = fullfile(paths.MNIparcs,'MNI_templates/MNI152_T1_1mm_brain.nii.gz');
        fileIn = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    else
        fileRef = fullfile(paths.MNIparcs,'MNI_templates/MNI152_T1_1mm.nii.gz');
        fileIn = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
    end
    if exist(fileIn,'file') && exist(fileRef,'file')
        fileMat_dof6 = fullfile(paths.T1.reg,'T12MNI_dof6.mat');
        fileOut = fullfile(paths.T1.reg,'T1_dof6');
        disp('          dof6')
        % Linear rigid body registration T1 to MNI
        sentence = sprintf('%s/flirt -ref %s -in %s -omat %s -out %s -cost %s -dof 6 -interp spline',...
            paths.FSL,fileRef,fileIn,fileMat_dof6,fileOut,configs.T1.flirtdof6cost);
         [~,result] = system(sentence); %#ok<*ASGLU>

    % inverse matrix flirt dof 6
        if exist(fileMat_dof6,'file')
            fileMat_dof6_inv = fullfile(paths.T1.reg,'MNI2T1_dof6.mat');
            sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',...
                paths.FSL,fileMat_dof6_inv,fileMat_dof6);
            [~,result] = system(sentence);
            if exist(fileMat_dof6_inv,'file') ~= 2
                warning('%s not generated. Exiting...',fileMat_dof6_inv)
                return
            end
        else
            warning('%s not generated. Exiting...',fileMat_dof6)
            return
        end
    else
        warning('%s and/or %s not found',fileIn,fileRef)
        return
    end

%% flirt dof 12
    fileIn = fullfile(paths.T1.reg,'T1_dof6.nii.gz');
    if exist(fileIn,'file') && exist(fileRef,'file')
        fileMat_dof12 = fullfile(paths.T1.reg,'T12MNI_dof12.mat');
        fileOut = fullfile(paths.T1.reg,'T1_dof12');
        disp('          dof12')
        % Linear affnie registration of T1 to MNI
        sentence = sprintf('%s/flirt -ref %s -in %s -omat %s -out %s -dof 12 -interp spline',...
            paths.FSL,fileRef,fileIn,fileMat_dof12,fileOut);%
        [~,result] = system(sentence);

        niifileOut = strcat(fileOut,'.nii.gz');
        if exist(fileMat_dof12,'file') ~= 2 || exist(niifileOut,'file') ~= 2
            warning('%s or %s not created. Exiting...',fileMat_dof12,niifileOut)
            return
        end

        % inverse matrix flirt dof 12
        fileMat_dof12_inv = fullfile(paths.T1.reg,'MNI2T1_dof12.mat');
        sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',...
            paths.FSL,fileMat_dof12_inv,fileMat_dof12);
        [~,result] = system(sentence);

        if exist(fileMat_dof12_inv,'file') ~= 2 
            warning('%s not created. Exiting...',fileMat_dof12_inv)
            return
        end  
    end    
%% fnirt
    fileIn = fullfile(paths.T1.reg,'T1_dof12.nii.gz');

    if exist(fileIn,'file') ~=2 || exist(fileRef,'file') ~= 2
        warning('%s or %s not found. Exiting...',fileIn,fileRef)
        return
    end

    fileOut = fullfile(paths.T1.reg,'T1_warped');
    fileWarp = fullfile(paths.T1.reg,'T12MNI_warp');
    disp('          nonlinear')
    % Nonlinear warp of T1 to MNI
    sentence = sprintf('%s/fnirt --ref=%s --in=%s --cout=%s --iout=%s --subsamp=%s',...
        paths.FSL,fileRef,fileIn,fileWarp,fileOut,configs.T1.fnirtSubSamp);
    [~,result] = system(sentence);
    disp(result)  
    niifileOut = strcat(fileOut,'.nii.gz');
    niifileWarp = strcat(fileWarp,'.nii.gz');
    if exist(niifileOut,'file') ~= 2 || exist(niifileWarp,'file') ~= 2
        warning('%s and/or %s not created. Exiting...',fileOut,fileWarp)
        return
    end

    % inverse warp fnirt
    fileRef = fullfile(paths.T1.reg,'T1_dof12');
    fileWarpInv = fullfile(paths.T1.reg,'MNI2T1_warp.nii.gz');
    sentence = sprintf('%s/invwarp --ref=%s --warp=%s --out=%s',...
        paths.FSL,fileRef,fileWarp,fileWarpInv);
    [~,result] = system(sentence);

    if exist(fileWarpInv,'file') ~= 2 
        warning('%s not created. Exiting...',fileWarpInv)
        return
    end
end

%% Transform parcellations from MNI to native subject space
% for every parcellation named in the batch setup +1 (Ventricle mask)
    numParcs=length(parcs.pdir); % EJC 2018.03.21 plabel replaced with pdir to avoid error when batching multiple subjects.
for k = 1:numParcs+1
    if numParcs+1 == k
        fileIn=(fullfile(paths.MNIparcs,'MNI_templates','MNI152_T1_1mm_VentricleMask.nii.gz'));
        parcs.plabel(k).name='CSFvent';
    elseif isempty(parcs.plabel)
        warning('No parcellations have been specified in batch setup')
        disp('....................................')
        warning('Proceeding with default shen 278 region parcellation')
        parcs.label(1).name='shen_org';
        parcs.pdir(1).name='shen_MNI152_org';
        fileIn=fullfile(paths.MNIparcs,parcs.pdir(1).name,strcat(parcs.pdir(1).name,'.nii.gz'));
    else
        fileIn=fullfile(paths.MNIparcs,parcs.pdir(k).name,strcat(parcs.pdir(k).name,'.nii.gz'));
    end
    
    if exist(fileIn,'file') ~=2
        warning('%s does not exist. Exiting..',fileIn)
        return
    end

    sprintf('     %s --> T1',parcs.plabel(k).name)
    disp('        unwarp')
    fileRef = fullfile(paths.T1.reg,'T1_dof12.nii.gz');
    fileOut = fullfile(paths.T1.reg,strcat(parcs.plabel(k).name,'_unwarped.nii.gz'));
    sentence = sprintf('%s/applywarp --ref=%s --in=%s --warp=%s --out=%s --interp=nn',...
        paths.FSL,fileRef,fileIn,fileWarpInv,fileOut);
    [~,result] = system(sentence);
    
    if exist(fileOut,'file') ~= 2
        warning('%s not created. Exiting..',fileOut)
        disp(result)
        return
    end
    
    % inv dof 12
    fileIn = fullfile(paths.T1.reg,strcat(parcs.plabel(k).name,'_unwarped.nii.gz'));
    fileRef = fullfile(paths.T1.reg,'T1_dof6.nii.gz');
    fileOut = fullfile(paths.T1.reg,strcat(parcs.plabel(k).name,'_unwarped_dof12.nii.gz'));
    disp('        dof12')
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
        paths.FSL,fileIn,fileRef,fileOut,fileMat_dof12_inv);
    [~,result] = system(sentence);
    
    if exist(fileOut,'file') ~= 2
        warning('%s not created. Exiting..',fileOut)
        disp(result)
        return
    end
    
    % inv dof 6
    fileIn = fullfile(paths.T1.reg,strcat(parcs.plabel(k).name,'_unwarped_dof12'));
    if configs.T1.useMNIbrain == 1
        fileRef = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    else
        fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
    end    
    fileOut = fullfile(paths.T1.reg,strcat(parcs.plabel(k).name,'_unwarped_dof12_dof6.nii.gz'));
    disp('        dof6')
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
        paths.FSL,fileIn,fileRef,fileOut,fileMat_dof6_inv);
    [~,result] = system(sentence);
    
    if exist(fileOut,'file') ~= 2
        warning('%s not created. Exiting..',fileOut)
        disp(result)
        return
    end
    
    if numParcs+1 == k 
        fileOut2 = fullfile(paths.T1.dir,'T1_mask_CSFvent.nii.gz');
    else
        fileOut2 = fullfile(paths.T1.dir,strcat('T1_parc_',parcs.plabel(k).name,'.nii.gz'));
        
    end
    sentence = sprintf('cp %s %s',fileOut,fileOut2);
    [status,result] = system(sentence);
    
    if status ~= 0
        warning('%s not created. Exiting..',fileOut2)
        disp(result)
        return
    end        
end
