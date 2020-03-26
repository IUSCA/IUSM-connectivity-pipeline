function [paths,flags,configs]=f_T1_prepare_A(paths,flags,configs)
%                            F_T1_PREPARE_A
% T1 precprocessing code that results in denoised/field-of-view cropped T1
% image, generated from provided dicom data. A brain mask as well as a
% brain extracted image are also created.
%
% This function is executed as part of the connectivity pipeline and is not
% mean to run stand-alone.
%
% Contributors:
%   Joaquin Goni, Purdue University
%   Joey Contreras, University Southern California
%   Mario Dzemidzic, Indiana University School of Medicine
%   Evgeny Chumin, Indiana University School of Medicine
%

%% Dicom to NiFTI
if flags.T1.dcm2niix==1
    % set dicom directy and identify file extension
    paths.T1.dcm = fullfile(paths.T1.dir,configs.name.dcmFolder);
    [dcm_ext]=find_dcm_ext(paths.T1.dcm);
    if isempty(dcm_ext)
        warning('No dicom (.IMA or .dcm) images found. Skipping further analysis')
        return
    end
    list_nii = dir(fullfile(paths.T1.dir,'*.nii*'));
    % Remove existing nifti images.
    if size(list_nii,1) >= 1
        sentence = sprintf('rm %s/*.nii*',paths.T1.dir);
        [~,result] = system(sentence); %#ok<*ASGLU>
    end
    % Convert dicom to nifti.
    disp('Converting Dicom-to-nifti')
    fileLog = sprintf('%s/dcm2niix.log',paths.T1.dir);
    sentence = sprintf('%s/connectome_scripts/dcm2niix/dcm2niix -f %s -o %s -v y -x y -b y %s > %s',paths.scripts,configs.name.T1,paths.T1.dir,paths.T1.dcm,fileLog);
    [status,result] = system(sentence);
    if status == 0
        sentence = sprintf('mv %s/%s.nii %s/%s_orig.nii',paths.T1.dir,configs.name.T1,paths.T1.dir,configs.name.T1);
        [~,result] = system(sentence); % Rename T1 to T1_orig
        if configs.T1.useCropped == 1
            sentence = sprintf('gzip -f %s/%s_Crop_1.nii',paths.T1.dir,configs.name.T1);
            [~,result] = system(sentence); % Gzip T1_Crop_1
            sentence = sprintf('mv %s/%s_Crop_1.nii.gz %s/%s.nii.gz',paths.T1.dir,configs.name.T1,paths.T1.dir,configs.name.T1);
            [~,result] = system(sentence); % Rename T1_Crop_1 to T1
        else
            sentence = sprintf('gzip -f %s/%s_orig.nii',paths.T1.dir,configs.name.T1);
            [~,result]=system(sentence);
            sentence = sprintf('mv %s/%s_orig.nii.gz %s/%s.nii.gz',paths.T1.dir,configs.name.T1,paths.T1.dir,configs.name.T1);
            [~,result]=system(sentence);
        end
    else
        warning('Conversion failed.')
        return
    end
end

%% T1 denoiser
        % ONLM denoise T1 data in a fixed intensity range.
if flags.T1.denoiser==1
    fileIn = fullfile(paths.T1.dir,sprintf('%s.nii.gz',configs.name.T1));
    if exist(fileIn,'file')
        disp('Starting T1 denoising')
        T1_fov = MRIread(fileIn);
        vol_denoised = f_denoise_T1(T1_fov.vol);
        T1_fov.vol = vol_denoised;
        MRIwrite(T1_fov,fullfile(paths.T1.dir,'T1_denoised.nii'));
    else
        warning(' %s.nii.gz not found. Exiting...',configs.name.T1)
        return
    end
end
if configs.T1.skipDenoise == 1
    fileIn = fullfile(paths.T1.dir,sprintf('%s.nii.gz',configs.name.T1));
    fileOut = fullfile(paths.T1.dir,'T1_denoised.nii.gz');
    [~,result]=system(sprintf('cp %s %s',fileIn,fileOut));
    disp(result)
    [~,result]=system(sprintf('gunzip -f %s',fileOut));
    disp(result)
    disp('Skipped Denoising. Copied T1 to T1_denoised for further processing.')
end

%% FSL ANAT
if flags.T1.anat==1
    paths.T1.anat = fullfile(paths.T1.dir,'T1_denoised.anat');
    disp('Running FSL_ANAT')
    % Remove anat directory if one already exists.
    if exist(paths.T1.anat,'dir')
        sentence = sprintf('rm -r %s',paths.T1.anat);
        [~,result] = system(sentence);
    end
    % Run FSL anat, for bias field correction and subcortical segmentation.
    fileIn1 = fullfile(paths.T1.dir,'T1_denoised.nii');
    if exist(fileIn1,'file')
        % strongbias should be more appropriate for multi-channel coils on 3T scanners.
        % add nocrop option if registration fails
        if configs.T1.bias == 0
            argbias = '--nobias';
        elseif configs.T1.bias == 1
            argbias = '--weakbias';
        else
            argbias = '';
        end
        if configs.T1.crop == 0
            argcrop = '--nocrop';
        else
            argcrop = '';
        end   
        arganat = sprintf('%s %s',argbias,argcrop);
        
        if configs.T1.padfix == 1
            FileAFNI=fullfile(paths.T1.dir,'pad5+orig.HEAD');
            sentence=sprintf('%s/3dZeropad -I 5 -prefix %s/pad5 %s; %s/3dAFNItoNIFTI -prefix %s/pad5 %s',...
                paths.AFNI,paths.T1.dir,fileIn1,paths.AFNI,paths.T1.dir,FileAFNI);
            [~,result]=system(sentence);
            fileIn1 = fullfile(paths.T1.dir,'pad5.nii');
            if exist(fileIn1,'file')
                [~,result]=system(sprintf('rm %s/pad5+orig*',paths.T1.dir));
            else
                fprintf(2,'zero padding T1_fov_denoised failed. Please debug..\n')
                return
            end
        end 
        sentence = sprintf('%s/fsl_anat --noreg --nononlinreg --noseg %s -i %s',paths.FSL,arganat,fileIn1);
        [~,result] = system(sentence);
        fileIn2 = fullfile(paths.T1.anat,'T1_biascorr.nii.gz');
        fileOut = fullfile(paths.T1.dir,'T1_fov_denoised.nii.gz');
        if exist(fileIn2,'file')
            % Unzip and copy the biascorrected image to the main derectory.
            sentence = sprintf('cp %s %s',fileIn2,fileOut);
            [status,result] = system(sentence);
            if status == 0
                sentence = sprintf('gunzip -f %s',fileOut);
                [status,result] = system(sentence);
                if status == 0
                    disp('Completed FSL_ANAT')
                else
                    warning('Gunzip failed. Check file ownership.') 
                    return
                end
            else
                warning('%s not created. Exiting...',fileOut)
                return
            end
        else
            warning('%s not found. Exiting...',fileIn2)
            return
        end
        % Subcortical masks
        disp('Copying subcortical segmentation.')
        fileIn = fullfile(paths.T1.anat,'T1_subcort_seg.nii.gz');
        fileOut = fullfile(paths.T1.dir,'T1_subcort_seg.nii.gz');
        if exist(fileIn,'file') ~= 2
            fprintf(2,'Subcortical segmentation not found. Exiting...\n')
            return
        end       
        sentence = sprintf('cp %s %s',fileIn,fileOut);
        [~,result] = system(sentence);
        
        % Generate a cerebellum mask using FSL's FIRST.
        if configs.T1.padfix == 1
            FileIn=fullfile(paths.T1.dir,'T1_fov_denoised.nii');
            FileAFNI=fullfile(paths.T1.dir,'pad5+orig.HEAD');
            sentence=sprintf('%s/3dZeropad -I 5 -prefix %s/pad5 %s; %s/3dAFNItoNIFTI -prefix %s/pad5 %s',...
                paths.AFNI,paths.T1.dir,FileIn,paths.AFNI,paths.T1.dir,FileAFNI);
            [~,result]=system(sentence);
            FileIn = fullfile(paths.T1.dir,'pad5.nii');
            if exist(FileIn,'file')
                [~,result]=system(sprintf('rm %s/pad5+orig*',paths.T1.dir));
            else
                fprintf(2,'zero padding T1_fov_denoised failed. Please debug..\n')
                return
            end
        else
            FileIn=fullfile(paths.T1.dir,'T1_fov_denoised.nii');
        end
        
        FileRoot=fullfile(paths.T1.dir,'subj_2_std_subc');
        FileMat=fullfile(paths.T1.dir,'subj_2_std_subc_cort.mat');
        FileOut1=fullfile(paths.T1.dir,'L_cerebellum.nii.gz');
        FileOut2=fullfile(paths.T1.dir,'R_cerebellum.nii.gz');
        paths.FSLroot=paths.FSL(1:end-4);
        FileModel1=fullfile(paths.FSLroot,'data/first/models_336_bin/intref_puta/L_Cereb.bmv');
        FileModel2=fullfile(paths.FSLroot,'data/first/models_336_bin/intref_puta/R_Cereb.bmv');
        FileRef1=fullfile(paths.FSLroot,'data/first/models_336_bin/05mm/L_Puta_05mm.bmv');
        FileRef2=fullfile(paths.FSLroot,'data/first/models_336_bin/05mm/R_Puta_05mm.bmv');
        sentence=sprintf('%s/first_flirt %s %s -cort',paths.FSL,FileIn,FileRoot);
        [~,result]=system(sentence);
        sentence=sprintf('%s/run_first -i %s -t %s -o %s -n 40 -m %s -intref %s',paths.FSL,FileIn,FileMat,FileOut1,FileModel1,FileRef1);
        [~,result]=system(sentence);
        sentence=sprintf('%s/run_first -i %s -t %s -o %s -n 40 -m %s -intref %s',paths.FSL,FileIn,FileMat,FileOut2,FileModel2,FileRef2);
        [~,result]=system(sentence);
        % Clean up the edges of the cerebellar mask.
        sentence=sprintf('%s/first_boundary_corr -s %s -i %s -b fast -o %s',paths.FSL,FileOut1,FileIn,FileOut1);
        [~,result]=system(sentence);
        sentence=sprintf('%s/first_boundary_corr -s %s -i %s -b fast -o %s',paths.FSL,FileOut2,FileIn,FileOut2);
        [~,result]=system(sentence);
        % Add the left and right cerebellum masks together.
        FileOut=fullfile(paths.T1.dir,'Cerebellum_bin.nii.gz');
        sentence=sprintf('%s/fslmaths %s -add %s %s',paths.FSL,FileOut1,FileOut2,FileOut);
        [~,result]=system(sentence);
        % remove extra slices if nesessary
        if configs.T1.padfix == 1
            FileAFNI=fullfile(paths.T1.dir,'cut5+orig.HEAD');
            sentence=sprintf('%s/3dZeropad -I -5 -prefix %s/cut5 %s; %s/3dAFNItoNIFTI -prefix %s/cut5 %s',...
                paths.AFNI,paths.T1.dir,FileOut,paths.AFNI,paths.T1.dir,FileAFNI);
            [~,result]=system(sentence);
            FileOut2 = fullfile(paths.T1.dir,'cut5.nii');
            if exist(FileOut2,'file')
                [~,result]=system(sprintf('rm %s/cut5+orig*',paths.T1.dir));
            else
                fprintf(2,'zero padded slice removal in Cerebellum_bin failed. Please debug..\n')
                return
            end
            [~,result]=system(sprintf('mv %s %s',FileOut2,FileOut));
            [~,result]=system(['rm ' FileOut2]);
        end 
        %-----------------------------------------------------------------%
        % Fill holes in the mask.
        sentence=sprintf('%s/fslmaths %s -fillh %s',paths.FSL,FileOut,FileOut);
        [~,result]=system(sentence);
        % Invert the cerebellum mask.
        FileInv=fullfile(paths.T1.dir,'Cerebellum_Inv.nii.gz');
        sentence=sprintf('%s/fslmaths %s -binv %s',paths.FSL,FileOut,FileInv);
        [~,result]=system(sentence);
        %-----------------------------------------------------------------%
        % 07.25.2017 EJC Remove intermediates of the clean-up.
        sentence = sprintf('rm %s*;rm %s*;rm %s*;',FileRoot,FileOut1,FileOut2);
        [~,result]=system(sentence);
        sentence = sprintf('rm %s/L_cerebellum_*;rm %s/R_cerebellum_*;',paths.T1.dir,paths.T1.dir);
        [~,result]=system(sentence); 
        if configs.T1.padfix == 1
            [~,result]=system(sprintf('rm %s/pad5.nii',paths.T1.dir));
        end
    else
        warning('%s not found. Exiting...',fileIn1)
        return
    end
end

%% T1 bet
if flags.T1.bet==1
    disp('Brain Extraction and Masking')
    fileIn = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
%----------    fileOut = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    fileOutroot = fullfile(paths.T1.dir,'T1_');
    switch configs.T1.antsTemplate
        case 'MICCAI'
            fileTemplate = fullfile(paths.scripts,'connectome_scripts/templates/brainmask_templates/MICCAI2012-Multi-Atlas-Challenge-Data/T_template0.nii.gz');
            fileProbability = fullfile(paths.scripts,'connectome_scripts/templates/brainmask_templates/MICCAI2012-Multi-Atlas-Challenge-Data/T_template0_BrainCerebellumProbabilityMask.nii.gz');
            fprintf('%s brain mask template selected\n',configs.T1.antsTemplate)
        case 'NKI'
            fileTemplate = fullfile(paths.scripts,'connectome_scripts/templates/brainmask_templates/NKI/T_template.nii.gz');
            fileProbability = fullfile(paths.scripts,'connectome_scripts/templates/brainmask_templates/NKI/T_template_BrainCerebellumProbabilityMask.nii.gz');
            fprintf('%s brain mask template selected\n',configs.T1.antsTemplate)
        case 'bet'
            fprintf('Using bet -f and -g inputs to perform fsl bet with -B option\n')
        otherwise
            fprintf('Unknown brain mask template selection: %s. Exiting...\n',configs.T1.antsTemplate)
    end   
            
    if exist(fileIn, 'file')
        fileIn2 = fullfile(paths.T1.dir,'T1_brain_mask.nii.gz');
        fileOut = fullfile(paths.T1.dir,'T1_brain.nii.gz');
        if strcmp(configs.T1.antsTemplate,'bet') == 1
            sentence = sprintf('%s/bet %s %s -B -m -f %.4f -g %.4f',...
                paths.FSL,fileIn,fileOut,configs.T1.betF,configs.T1.betG);
            [status,result] = system(sentence);
        else
            % ANTS brain extraction
            ANTSlog = fullfile(paths.T1.dir,'ants_bet.log');
            sentence = sprintf('%s/antsBrainExtraction.sh -d 3 -a %s -e %s -m %s -o %s > %s',...
                paths.ANTS,fileIn,fileTemplate,fileProbability,fileOutroot,ANTSlog);
            [status,result]=system(sentence);
            [status,result]=system(sprintf('mv %s/T1_BrainExtractionMask.nii.gz %s',paths.T1.dir,fileIn2));
            [status,result]=system(sprintf('mv %s/T1_BrainExtractionBrain.nii.gz %s',paths.T1.dir,fileOut));
        end
        fileOut2 = fullfile(paths.T1.dir,'T1_brain_mask_filled.nii.gz');
        if status == 0 && exist(fileIn2, 'file')
            % Fill holes in the brain mask.
            sentence = sprintf('%s/fslmaths %s -fillh %s',paths.FSL,fileIn2,fileOut2);
            [status,result] = system(sentence);
            if status == 0 && exist(fileOut2, 'file')
                disp('Bet completed.')
            else
                warning('%s not created. Exiting...',fileOut2)
                return
            end
        else
            warning('%s not found. Exiting...',fileIn2)
            return
        end
    else
        warning('%s not found. Exiting...',fileIn)
        return
    end
end

%% T1 Brain Re-Extract
  % Use the filled brain mask to extract the brain.
if flags.T1.re_extract ==1
    fileIn = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
    fileMask = fullfile(paths.T1.dir,'T1_brain_mask_filled.nii.gz');
    fileOut = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    if exist(fileIn,'file') && exist(fileMask,'file')
        sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn,fileMask,fileOut);
        [status,result] = system(sentence);
        if status == 0 && exist(fileOut,'file')
            strdisp = strcat(fileOut,' created.');
            disp(strdisp)
        end
    else
        warning('%s and/or %s not found. Exiting...',fileIn,fileMask)
        return
    end
end

