function [paths,flags,configs,parcs]=f_T1_prepare_B(paths,flags,configs,parcs)
%                               F_T1_PREPARE_B
% Registration of subject anatomical T1 image to MNI with subsequent
% transformation of parcellation images into subject native space.
%
% Segmetation of subject T1 into tissue-types.
%
% Contributors:
%   Joaquin Goni, Purdue University
%   Joey Contreras, University Southern California
%   Mario Dzemidzic, Indiana University School of Medicine
%   Evgeny Chumin, Indiana University School of Medicine

%% Registration of subjec to MNI
if flags.T1.reg2MNI==1
    % Transform subject T1 into MNI space, then using inverse matrices
    % transform yeo7, yeo17, shen parcellations and MNI ventricles mask
    % into subject native space.
    disp('Registration between Native T1 and MNI space') 
    [paths,configs,parcs]=fsl_registration_parcellations_non_linear(paths,configs,parcs);
end
paths.MNIparcs = fullfile(paths.scripts,'connectome_scripts/templates/MNIparcs');
%% Tissue-type segmentation; cleaning; and gray matter masking of parcellations
if flags.T1.seg==1
    disp('Tissue-type Segmentation') 
    % Check that T1_brain image indeed exists.
    fileIn = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    if exist(fileIn,'file') ~=2 
        warning('%s not found. Exiting...',fileIn')
        return
    end
    
    % FSL fast tissue-type segmentation (GM, WM, CSF)
    sentence = sprintf('%s/fast -H %s %s',paths.FSL,configs.T1.segfastH,fileIn);
    [~,result] = system(sentence); %#ok<*ASGLU>

    %% CSF masks
    fileIn = fullfile(paths.T1.dir,'T1_brain_seg.nii.gz');
    fileOut = fullfile(paths.T1.dir,'T1_CSF_mask');
    if exist(sprintf(fileIn),'file') ~= 2
        warning('%s not found. Exiting...',fileIn')
        return
    end
        % one is equal to CSF
    sentence = sprintf('%s/fslmaths %s -thr 1 -uthr 1 %s',...
        paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    sentence = sprintf('%s/fslmaths %s/T1_CSF_mask.nii.gz -mul -1 -add 1 %s/T1_CSF_mask_inv.nii.gz',...
        paths.FSL,paths.T1.dir,paths.T1.dir);
    [~,result] = system(sentence);  
    
    %% Subcortical masks
    fileSub = fullfile(paths.T1.dir,'T1_subcort_seg.nii.gz');
    if exist(fileSub,'file') ~= 2
        warning('%s not found. Exiting...',fileSub)
        return
    end
    
    %binarize subcortical segmentation to a mask
    fileIn = fullfile(paths.T1.dir,'T1_subcort_mask.nii.gz');
    sentence = sprintf('%s/fslmaths %s -bin %s',paths.FSL,fileSub,fileIn);
    [~,result] = system(sentence);
    % REmove CSF contamination
    fileMas = fullfile(paths.T1.dir,'T1_CSF_mask_inv.nii.gz');
    sentence = sprintf('%s/fslmaths %s -mas %s %s',...
        paths.FSL,fileIn,fileMas,fileIn);
    [~,result] = system(sentence);
    sentence = sprintf('%s/fslmaths %s -mas %s %s',...
        paths.FSL,fileSub,fileMas,fileSub);
    [~,result] = system(sentence);

    sentence = sprintf('%s/fslmaths %s -mul -1 -add 1 %s/T1_subcort_mask_inv.nii.gz',...
        paths.FSL,fileIn,paths.T1.dir);
    [~,result] = system(sentence);
    
    %% Adding FIRST subcortical into tissue segmentation
    sentence = sprintf('%s/fslmaths %s/T1_brain_seg -mul %s/T1_subcort_mask_inv %s/T1_brain_seg_best',...
        paths.FSL,paths.T1.dir,paths.T1.dir,paths.T1.dir);
    [~,result] = system(sentence);

    sentence = sprintf('%s/fslmaths %s/T1_subcort_mask -mul 2 %s/T1_subcort_seg_add',...
        paths.FSL,paths.T1.dir,paths.T1.dir);
    [~,result] = system(sentence);

    sentence = sprintf('%s/fslmaths %s/T1_brain_seg_best -add %s/T1_subcort_seg_add %s/T1_brain_seg_best',...
        paths.FSL,paths.T1.dir,paths.T1.dir,paths.T1.dir);
    [~,result] = system(sentence);

    %% Separating Tissue types
    listTissue = {'CSF','GM','WM'};
    fileIn = fullfile(paths.T1.dir,'T1_brain_seg_best.nii.gz');
    for i=1:3 % 1=CSF 2=GM 3=WM
        fileOut = fullfile(paths.T1.dir,sprintf('T1_%s_mask',listTissue{i}));
        sentence = sprintf('%s/fslmaths %s -thr %d -uthr %d -div %d %s',...
            paths.FSL,fileIn,i,i,i,fileOut);
        [~,result] = system(sentence);
        %Erode each tissue mask
        sentence=sprintf('%s/fslmaths %s -ero %s_eroded',paths.FSL,fileOut,fileOut);
        [~,result]=system(sentence);
        if i == 3 % if WM
            WMeroded = fullfile(paths.T1.dir,'T1_WM_mask_eroded.nii.gz');
            % 2nd WM erosion 
            sentence = sprintf('%s/fslmaths %s -ero %s',paths.FSL,WMeroded,WMeroded);
            [~,result]=system(sentence);
            % 3rd WM erosion
            sentence = sprintf('%s/fslmaths %s -ero %s',paths.FSL,WMeroded,WMeroded);
            [~,result]=system(sentence);
        end
    end
    
    % apply as CSF ventricles mask 
    fileIn = fullfile(paths.T1.dir,'T1_CSF_mask_eroded.nii.gz');
    fileOut = fullfile(paths.T1.dir,'T1_CSFvent_mask_eroded');
    fileMas = fullfile(paths.T1.dir,'T1_mask_CSFvent.nii.gz');
    if exist(fileMas,'file') ~= 2
        warning('%s not found. Exiting...',fileMas')
        return
    end
    sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileIn,fileMas,fileOut);
    [~,result] = system(sentence);
    
%% WM CSF sandwich
     disp('WM/CSF sandwich')
    % Remove any gray matter voxels that are withing one dilation of CSF and white matter.
    % Dilate WM mask
     fileIn = fullfile(paths.T1.dir,'T1_WM_mask.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_WM_mask_dil');
     sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileIn,fileOut);
     [~,result] = system(sentence);

     % Dilate CSF mask
     fileIn = fullfile(paths.T1.dir,'T1_CSF_mask.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_CSF_mask_dil');
     sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileIn,fileOut);
     [~,result] = system(sentence);
     
     % Add the dilated masks together.
     fileIn1 = fullfile(paths.T1.dir,'T1_WM_mask_dil.nii.gz');
     fileIn2 = fullfile(paths.T1.dir,'T1_CSF_mask_dil.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     sentence = sprintf('%s/fslmaths %s -add %s %s',paths.FSL,fileIn1,fileIn2,fileOut);
     [~,result] = system(sentence);
     
     % Theshold the image at 2, isolating WM, CSF interface.
     fileIn = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     sentence = sprintf('%s/fslmaths %s -thr 2 %s',paths.FSL,fileIn,fileOut);
     [~,result] = system(sentence);
     
     % Multiply the interface by the native space ventricle mask.
     fileIn1 = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     fileIn2 = fullfile(paths.T1.dir,'T1_mask_CSFvent.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn1,fileIn2,fileOut);
     [~,result] = system(sentence);
     
     % Using fsl cluster identify the largest contiguous cluster, and save
     % it out as a new mask.
     filein = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     textout = fullfile(paths.T1.dir,'WM_CSF_sandwich_clusters.txt');
     sentence = sprintf('%s/cluster --in=%s --thresh=1 --osize=%s >%s',paths.FSL,filein,filein,textout);
     [~,result]=system(sentence);
     cluster=dlmread(fullfile(paths.T1.dir,'WM_CSF_sandwich_clusters.txt'),'\t',[1 1 1 1]);
     sentence = sprintf('%s/fslmaths %s -thr %d %s',paths.FSL,filein,cluster,filein);
     [~,result]=system(sentence);
     
     % Binarize and invert the single cluster mask.
     fileIn = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     fileOut = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich');
     sentence = sprintf('%s/fslmaths %s -binv %s',paths.FSL,fileIn,fileOut);
     [~,result] = system(sentence);
     
     % Filter the GM mask with obtained CSF_WM sandwich.
     fileIn1 = fullfile(paths.T1.dir,'T1_WM_CSF_sandwich.nii.gz');
     fileIn2 = fullfile(paths.T1.dir,'T1_GM_mask.nii.gz');
     fileOut = fileIn2;
     sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn1,fileIn2,fileOut);
     [~,result] = system(sentence);
end

%% Intersect of parcellations with GM
if flags.T1.parc==1
% Gray matter masking of native space parcellations
    counter = 0;
    for k=1:length(parcs.pdir) % for every indicated parcellation
        fileIn = fullfile(paths.T1.dir,strcat('T1_parc_',parcs.plabel(k).name,'.nii.gz'));
        if exist(fileIn,'file') ~= 2
            warning('%s not found. Exiting...',fileIn')
            return
        end
        fprintf('%s parcellation intersection with GM\n',parcs.plabel(k).name) 
        fileOut = fullfile(paths.T1.dir,strcat('T1_parc_',parcs.plabel(k).name,'_dil.nii.gz'));
        % Dilate the parcellation.
        sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileIn,fileOut);
        [~,result] = system(sentence);
    
        % Iteratively mask the dilated parcellation with GM.
        fileMul = fullfile(paths.T1.dir,'T1_GM_mask.nii.gz');
        if exist(fileMul,'file') ~= 2
            warning('%s not found. Exiting...',fileMul')
            return 
        end
        % Apply subject GM mask
        fileOut2 = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
        sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileOut,fileMul,fileOut2);
        [~,result] = system(sentence);
        % Dilate and remask to fill GM mask a set number of times
        fileOut3 = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'_dil.nii.gz'));
        for i=1:configs.T1.numDilReMask
            sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileOut2,fileOut3);
            [~,result]=system(sentence);
            sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileOut3,fileMul,fileOut2);
            [~,result]=system(sentence);
        end
        % 07.25.2017 EJC Remove the left over dil parcellation images.
        sentence = sprintf('rm %s %s %s',fileOut,fileOut3);
        [~,result]=system(sentence);
        if parcs.pcort(k).true == 1
            counter=counter+1;
            %-------------------------------------------------------------------------%
            % Clean up the cortical parcellation by removing subcortical and
            % cerebellar gray matter.
            paths.T1.reg=fullfile(paths.T1.dir,'registration');
            if counter == 1
                % Generate inverse subcortical mask to isolate cortical portion of parcellation.
                fileIn = fullfile(paths.T1.dir,'T1_subcort_mask.nii.gz');
                fileOut = fullfile(paths.T1.dir,'T1_subcort_mask_dil.nii.gz');
                fileMas = fullfile(paths.T1.dir,'T1_GM_mask.nii.gz');
                sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileIn,fileOut);
                [~,result]=system(sentence);
                sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileOut,fileMas,fileOut);
                [~,result]=system(sentence);
                fileMas2 = fullfile(paths.T1.dir,'T1_subcort_mask_dil_inv.nii.gz');
                sentence = sprintf('%s/fslmaths %s -binv %s',paths.FSL,fileOut,fileMas2);
                [~,result]=system(sentence);
            end
            %---------------------------------------------------------%
            % Apply subcortical inverse to cortical parcellations.
            fileOut = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileOut,fileMas2,fileOut);
            [~,result]=system(sentence);
            %---------------------------------------------------------%
            % Generate a cerebellum mask using FSL's FIRST.
            if counter == 1
                % inverse transfrom the MNI cerebellum mask
                fileIn = fullfile(paths.MNIparcs,'MNI_templates/MNI152_T1_cerebellum.nii.gz');
                fileWarpInv = fullfile(paths.T1.reg,'MNI2T1_warp.nii.gz');
                fileRef = fullfile(paths.T1.reg,'T1_dof12.nii.gz');
                fileOut = fullfile(paths.T1.reg,'cerebellum_unwarped.nii.gz');
                sentence = sprintf('%s/applywarp --ref=%s --in=%s --warp=%s --out=%s --interp=nn',...
                    paths.FSL,fileRef,fileIn,fileWarpInv,fileOut);
                [~,result] = system(sentence);

                fileIn = fullfile(paths.T1.reg,'cerebellum_unwarped.nii.gz');
                fileRef = fullfile(paths.T1.reg,'T1_dof6.nii.gz');
                fileOut = fullfile(paths.T1.reg,'cerebellum_unwarped_dof12.nii.gz');
                fileMat_dof12_inv = fullfile(paths.T1.reg,'MNI2T1_dof12.mat');
                sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
                    paths.FSL,fileIn,fileRef,fileOut,fileMat_dof12_inv);
                [~,result] = system(sentence);

                fileMat_dof6_inv = fullfile(paths.T1.reg,'MNI2T1_dof6.mat');
                fileIn = fullfile(paths.T1.reg,'cerebellum_unwarped_dof12');
                if configs.T1.useMNIbrain == 1
                    fileRef = fullfile(paths.T1.dir,'T1_brain.nii.gz');
                else
                    fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
                end    
                fileOut = fullfile(paths.T1.reg,'Cerebellum_bin.nii.gz');
                sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
                    paths.FSL,fileIn,fileRef,fileOut,fileMat_dof6_inv);
                [~,result] = system(sentence);

                % Generate a cerebellar inverse mask
                fileIn = fullfile(paths.T1.reg,'Cerebellum_bin.nii.gz');
                fileOut = fullfile(paths.T1.dir,'Cerebellum_Inv.nii.gz');
                [~,result]=system(sprintf('%s/fslmaths %s -binv %s',paths.FSL,fileIn,fileOut));
            end
        %-------------------------------------------------------------------------%    
        % Remove any parcellation contamination of the cerebellum.
        FileIn=fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
        FileInv=fullfile(paths.T1.dir,'Cerebellum_Inv.nii.gz');
        sentence=sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,FileIn,FileInv,FileIn);
        [~,result]=system(sentence);
        %-------------------------------------------------------------------------%  
        end
        %% add subcortical fsl parcellation to cortical parcellations
        if configs.T1.addsubcort == 1
            fileSubcort = fullfile(paths.T1.dir,'T1_subcort_seg.nii.gz');
            volParc=MRIread(FileIn);
            MaxID = max(max(max(volParc.vol)));
            volSubcort=MRIread(fileSubcort);
            volSubcort.vol(volSubcort.vol==16)=0;
            %----------------------------------------------------------%
            if parcs.pnodal(k).true == 1 
                ids=unique(volSubcort.vol);
                for s=1:length(ids)
                    if ids(s)>0
                        volSubcort.vol(volSubcort.vol==ids(s))=MaxID+(s-1);
                    end
                end
            elseif parcs.pnodal(k).true == 0
                volSubcort.vol(volSubcort.vol>0)=MaxID+1;
            end
            %----------------------------------------------------------%
            subcorMask=volSubcort.vol > 0;
            volParc.vol(subcorMask)=0;
            volParc.vol=volParc.vol+volSubcort.vol;
            MRIwrite(volParc,FileIn)
        end
%-------------------------------------------------------------------------%
        % 07.26.2017 EJC Dilate the final GM parcellations. 
        % NOTE: These will be used by f_functional_connectivity to bring parcellations into epi space.
        fileOut4 = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'_dil.nii.gz'));
        sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileOut2,fileOut4);
        [~,result]=system(sentence);
        if ~isempty(result)
            warning('Dilation of %s parcellation error! See return below for details.',parcs.plabel(k).name);
            disp(result)
        end
    end
end 