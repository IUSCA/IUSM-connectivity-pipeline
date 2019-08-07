function [paths,flags,configs,parcs] = f_prob_conn(paths,flags,configs,parcs,params)
%                       F_PROB_CONN
% Probabilistic tractography with FSL bedpostX and probtrackX.

% Contributors:
%   Evgeny Chumin, Indiana University School of Medicine
%
paths.DWI.EDDY=fullfile(paths.DWI.dir,'EDDY');
if exist(paths.DWI.EDDY,'dir')
    paths.DWI.DTIfit=fullfile(paths.DWI.dir,'DTIfit');
    if ~exist(paths.DWI.DTIfit,'dir')
        warning('Path to DTIfit directory does not exist. Exiting...')
        return
    end
else
    warning('Path to EDDY directory does not exist. Exiting...')
    return
end

%% Prepare bedpostX directory
if flags.DWI.bedpostXprep == 1
    disp('------------')
    disp('FSL BedpostX')
    disp('------------')    
    disp('---------------------')
    disp('     Setting up Input')
    disp('---------------------')
    % Extract subject identifier
    tmp=extractBetween(paths.DWI.dir,paths.data,'/DWI');
    tmp=extractAfter(tmp,'/'); subID=tmp{1}; clear tmp %#ok<*NASGU>
    
    paths.DWI.bpxIN = fullfile(paths.DWI.dir,subID);
    if ~exist(paths.DWI.bpxIN,'dir')
        mkdir(paths.DWI.bpxIN) % create directory for bedpost inputs
    end 
    
    % list of bedpostx inputs
    fileNii=fullfile(paths.DWI.EDDY,'eddy_output.nii.gz');
    filebvecs=fullfile(paths.DWI.EDDY,'eddy_output.eddy_rotated_bvecs');
    filebvals=fullfile(paths.DWI.DTIfit,'3_DWI.bval');
    fileb0Mask=fullfile(paths.DWI.DTIfit,'b0_1st_mask.nii.gz');
    % link input files into the directory
    sentence=sprintf('ln -s %s %s/data.nii.gz #',fileNii,paths.DWI.bpxIN);
        [~,result]=system(sentence);
    if isempty(result) || contains(result,'File exists')==1
        disp('Link to data created')
    else
        warning(result)
    end
    sentence=sprintf('ln -s %s %s/bvecs #',filebvecs,paths.DWI.bpxIN);
        [~,result]=system(sentence);
    if isempty(result) || contains(result,'File exists')==1
        disp('Link to bvecs created')
    else
        warning(result)
    end    
    sentence=sprintf('ln -s %s %s/bvals #',filebvals,paths.DWI.bpxIN);
        [~,result]=system(sentence);
    if isempty(result) || contains(result,'File exists')==1
        disp('Link to bvals created')
    else
        warning(result)
    end    
    sentence=sprintf('ln -s %s %s/nodif_brain_mask.nii.gz #',fileb0Mask,paths.DWI.bpxIN);
        [~,result]=system(sentence);
    if isempty(result) || contains(result,'File exists')==1
        disp('Link to nodif_brain_mask created')
    else
        warning(result)
    end      
    
    % check that all inputs are present in the directory
    sentence=sprintf('%s/bedpostx_datacheck %s',paths.FSL,paths.DWI.bpxIN);
    [~,result]=system(sentence);
    disp(result)
    if contains(result,'does not exist')==1 || contains(result,'No such file or directory')==1
        warning('BedpostX inputs missing. Check that the following exist and links have been made')
        disp(fileNii); disp(filebvecs); disp(filebvals); disp(fileb0Mask)
        disp('Terminating subject ...'); return
    end

    % Check if processing will be outsources and copy files accordingly
    if configs.DWI.outsourceBpX == 1
        if ~isempty(paths.DWI.groupDir)
            if ~exist(paths.DWI.groupDir,'dir')
                mkdir(paths.DWI.groupDir)
            end
        sentence=sprintf('cp -L -R -f %s %s',paths.DWI.bpxIN,paths.DWI.groupDir);
        [~,result]=system(sentence);
        else
            warning('No destination group directory specified in batch file')
        end
    end
end

%% Run bedpostx orientation distribution estimation
if flags.DWI.runBedpostX == 1
    disp('---------------------------------------')
    disp('     Estimation of Diffusion Parameters')
    disp('---------------------------------------')
    % Extract subject identifier
    tmp=extractBetween(paths.DWI.dir,paths.data,'/DWI');
    tmp=extractAfter(tmp,'/'); subID=tmp{1}; clear tmp %#ok<*NASGU>
    % Check that all directories and inputs are present
    paths.DWI.bpxIN = fullfile(paths.DWI.dir,subID);
    if ~exist(paths.DWI.bpxIN,'dir')
        warning('BedpostX input directory does not exist. Exiting ...'); return
    end
    sentence=sprintf('%s/bedpostx_datacheck %s',paths.FSL,paths.DWI.bpxIN);
    [~,result]=system(sentence);
    disp(result)
    if contains(result,'does not exist')==1 || contains(result,'No such file or directory')==1
        warning('BedpostX inputs missing. Check that the links have been made')
        disp('Terminating subject ...'); return
    end
    
    % Set up Parallel environment
    if configs.parallel == 1
        [rtrn, nproc] = system('nproc'); % Get the number of available processors
        num_proc2use = int2str(str2num(nproc) * configs.UsageCPU); %#ok<ST2NM> % Determine the number of proc to use
        setenv('FSLPARALLEL',num_proc2use);
    end
    
    % Run
    disp('Starting bedpostx. This will take a while. Come back tomorrow.')
    
    sentence=sprintf('%s/bedpostx %s --nf=%.0f --fudge=%.0f --bi=%.0f'...
        ,paths.FSL,paths.DWI.bpxIN,params.DWI.fibres,params.DWI.weight,params.DWI.burnin);
    [~,result]=system(sentence);
    disp(result)
    disp('BedpostX finished!')
end

%% Registration of T1 to DWI
if flags.DWI.reg2DWI == 1
    % register T1 to DWI space
    disp('---------------------')
    disp('Registering T1 to DWI')
    disp('---------------------')
    
    disp('rigid body dof6 to T1')
    fileIn=fullfile(paths.DWI.DTIfit,'3_DWI_FA.nii.gz');
    fileRef=fullfile(paths.T1.dir,'T1_brain.nii.gz');
    fileMat2T1=fullfile(paths.DWI.dir,'FA_2T1dof6.mat');
    fileDof6=fullfile(paths.DWI.dir,'FA_dof6.nii.gz');
    sentence=sprintf('%s/flirt -in %s -ref %s -omat %s -dof 6 -cost mutualinfo -interp spline -out %s',...
        paths.FSL,fileIn,fileRef,fileMat2T1,fileDof6);
    [~,result]=system(sentence); %#ok<*ASGLU>
    
    disp('dof6 inverse')
    if exist(fileMat2T1,'file')
        fileMat2FA=fullfile(paths.DWI.dir,'T1_2FAdof6.mat');
        sentence=sprintf('%s/convert_xfm -omat %s -inverse %s',...
            paths.FSL,fileMat2FA,fileMat2T1);
        [~,result]=system(sentence);
    else
        warning('FA_2_T1 transformation does not exist. Check FLIRT output. Exiting...')
        return
    end
    
    if configs.DWI.Warp == 1  
        disp('warp to FA')
        fileWout=fullfile(paths.DWI.dir,'T12FA_warped.nii.gz');
        filewarp=fullfile(paths.DWI.dir,'T12FA_warpfield.nii.gz');
        sentence=sprintf('%s/fnirt --in=%s --ref=%s --aff=%s --iout=%s --cout=%s --subsamp=8,4,2,2',...
            paths.FSL,fileRef,fileIn,fileMat2FA,fileWout,filewarp);
        [~,result]=system(sentence);
    end
    
        disp('apply to T1 & tissue-type volumes')
        fileT1out=fullfile(paths.DWI.dir,'rT1_brain.nii.gz');
        fileGM=fullfile(paths.T1.dir,'T1_GM_mask.nii.gz');
        fileGMout=fullfile(paths.DWI.dir,'rT1_GM_mask.nii.gz');   
    if configs.DWI.Warp == 1
        sentence=sprintf('%s/applywarp -i %s -o %s -r %s -w %s --interp=spline',...
            paths.FSL,fileRef,fileT1out,fileIn,filewarp);
        [~,result]=system(sentence);
        sentence=sprintf('%s/applywarp -i %s -o %s -r %s -w %s --interp=nn',...
            paths.FSL,fileGM,fileGMout,fileIn,filewarp);
        [~,result]=system(sentence);
    else
        sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp spline',...
            paths.FSL,fileRef,fileIn,fileMat2FA,fileT1out);
        [~,result]=system(sentence);
        sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
            paths.FSL,fileGM,fileIn,fileMat2FA,fileGMout);
        [~,result]=system(sentence);
    end
    
        disp('apply to node parcellations')
        for k=1:length(parcs.pnodal)
            if parcs.pnodal(k).true == 1
                fileParcT1=fullfile(paths.T1.dir,sprintf('T1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
                fileParcDWI=fullfile(paths.DWI.dir,sprintf('rT1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
                if configs.DWI.Warp == 1
                    sentence=sprintf('%s/applywarp -i %s -o %s -r %s -w %s --interp=nn',...
                        paths.FSL,fileParcT1,fileParcDWI,fileIn,filewarp);
                    [~,result]=system(sentence);
                else
                    sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                        paths.FSL,fileParcT1,fileIn,fileMat2FA,fileParcDWI);
                    [~,result]=system(sentence);
                end
            end
        end
end

%% Run probtrackx probabilistic tractography
if flags.DWI.probtrackX == 1
    disp('--------------')
    disp('FSL ProbtrackX')
    disp('--------------')
    % Extract subject identifier
    tmp=extractBetween(paths.DWI.dir,paths.data,'/DWI');
    tmp=extractAfter(tmp,'/'); subID=tmp{1}; clear tmp %#ok<*NASGU>

    if configs.DWI.importBpX == 1
        disp('-------------------------------')
        disp('     Importing BedpostX Results')
        disp('-------------------------------')
        if ~exist(paths.DWI.groupDir,'dir')
            warning('Specified group bedpost directory does not exist. Exiting ...');return
        end
        
        importDir=fullfile(paths.DWI.groupDir,strcat(subID,'.bedpostX'));
        if ~exist(importDir,'dir')
            warning('Import Error: Subject bedpostX output directory not found. Exiting...');return
        end

        paths.DWI.bpxOUT = fullfile(paths.DWI.dir,strcat(subID,'.bedpostX'));
        sentence=sprintf('cp -fR %s %s',importDir,paths.DWI.bpxOUT);
        [~,result]=system(sentence);
    end
        
    % Check that the output directory exists  
    paths.DWI.bpxOUT = fullfile(paths.DWI.dir,strcat(subID,'.bedpostX'));
    if ~exist(paths.DWI.bpxOUT,'dir')
        warning('BedpostX output directory does not exist. Exiting ...'); return
    end
    
    % create output directory
    paths.DWI.nodeTracks=fullfile(paths.DWI.dir,'seedTracts');
    if ~exist(paths.DWI.nodeTracks,'dir')
       mkdir(paths.DWI.nodeTracks)
    end
    
    % path to T1 registration parameters
    paths.T1.reg=fullfile(paths.T1.dir,'registration');
    % Parse seed inputs
    switch params.DWI.SeedType(1) %seed type
        case 1 % parcellation node
            disp('-------------------------------')
            disp('     Working on: Nodes as Seeds')
            disp('-------------------------------')
            if ~isempty(configs.DWI.SeedIdx) 
                if exist(paths.DWI.seedLabels,'file')
                    % read in label file
                    fid=fopen(paths.DWI.seedLabels);
                    slabels=textscan(fid,'%s');
                    fclose(fid);
                    if length(slabels{1})~=length(configs.DWI.SeedIdx)
                        disp('Number of labels does not match number of ROI.')
                        disp('Double check the label text file for proper size/format.')
                        disp('......')
                        disp('Proceeding without labels. ROI will be identified by integer values')
                        slabelON=0;
                    else
                        slabelON=1;
                    end
                else
                    disp('Seed label file does not exist. Proceeding without labels.')
                    slabelON=0;
                end
                for s=1:length(configs.DWI.SeedIdx) % for every seed in directory
                    % generate seed volume
                    for k=1:length(parcs.pnodal)
                        if parcs.pnodal(k).true == 1
                            fileIn=fullfile(paths.DWI.dir,sprintf('rT1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
                            if slabelON==1
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('seed_%s_%s.nii.gz',parcs.plabel(k).name,slabels{1}{s}));
                            elseif slabelON==0
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('seed_%s_%d.nii.gz',parcs.plabel(k).name,configs.DWI.SeedIdx(s)));
                            end
                            sentence=sprintf('%s/fslmaths %s -thr %d -uthr %d %s',...
                                paths.FSL,fileIn,configs.DWI.SeedIdx(s),configs.DWI.SeedIdx(s),fileOut);
                            [~,result]=system(sentence);
                            % Store names of seed files for probtrackX use
                            seedFiles{s,1}=fileOut; %#ok<*AGROW>
                        end
                    end
                end
            else
                warning('configs.DWI.SeedIdx is empty. Exiting...')
                return
            end  
        case 2 %MNI seed node
            disp('-------------------------------')
            disp('          Working on: MNI Seeds')
            disp('-------------------------------')
            if exist(paths.DWI.MNIseedImage,'file')
                % Prepare the seed images
                % Tranform to diffusion space
                %concatenate linear transformations form T1
                MNI2T1_linear=fullfile(paths.T1.reg,'MNI2T1_linear.mat');
                if ~exist(MNI2T1_linear,'file')
                    dof6_2T1=fullfile(paths.T1.reg,'MNI2T1_dof6.mat');
                    dof12_2T1=fullfile(paths.T1.reg,'MNI2T1_dof12.mat');
                    sentence=sprintf('%s/convert_xfm -omat %s -concat %s %s',...
                        paths.FSL,MNI2T1_linear,dof12_2T1,dof6_2T1);
                    [~,result]=system(sentence);
                end
                % apply MNI2T1 transformation
                MNI2T1_warp=fullfile(paths.T1.reg,'MNI2T1_warp.nii.gz');
                fileT1ref=fullfile(paths.T1.dir,'T1_fov_denoised.nii');
                [~,name,ext]=fileparts(paths.DWI.MNIseedImage);
                T1seedImage=fullfile(paths.DWI.nodeTracks,sprintf('T1_%s',[name ext]));
                sentence=sprintf('%s/applywarp -i %s -r %s -o %s -w %s --postmat=%s --interp=nn',...
                    paths.FSL,paths.DWI.MNIseedImage,fileT1ref,T1seedImage,MNI2T1_warp,MNI2T1_linear);
                [~,result]=system(sentence);
                % apply T1 to DWI transformation
                fileMat2FA=fullfile(paths.DWI.dir,'T1_2FAdof6.mat');
                filewarp=fullfile(paths.DWI.dir,'T12FA_warpfield.nii.gz');
                DWIseedImage=fullfile(paths.DWI.nodeTracks,sprintf('DWI_%s',[name ext]));
                fileRef=fullfile(paths.DWI.DTIfit,'3_DWI_FA.nii.gz');
                if exist(filewarp,'file')
                    sentence=sprintf('%s/applywarp -i %s -o %s -r %s -w %s --interp=nn',...
                        paths.FSL,T1seedImage,DWIseedImage,fileRef,filewarp);
                    [~,result]=system(sentence);
                else
                    sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                        paths.FSL,T1seedImage,fileRef,fileMat2FA,DWIseedImage);
                    [~,result]=system(sentence);
                end

                % mask & dilate image within gray matter
                GM_Mask=fullfile(paths.DWI.dir,'rT1_GM_mask.nii.gz');
                fileMasked=fullfile(paths.DWI.nodeTracks,sprintf('DWI_GM_%s',[name ext]));
                sentence=sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,DWIseedImage,GM_Mask,fileMasked);
                [~,result]=system(sentence);
%                 sentence=sprintf('%s/fslmaths %s -dilD -mas %s %s',paths.FSL,fileMasked,GM_Mask,fileMasked);
%                 [~,result]=system(sentence);

                % Read in and sort seed image
                seeds=MRIread(fileMasked);
                num4rthDIM=size(seeds.vol,4);
                if num4rthDIM==1 %3D
                    imgType=3;
                    seedIds=unique(seeds.vol);
                    seedIds(seedIds==0)=[];
                    numROI=length(seedIds);   
                elseif num4rthDIM > 1 %4D
                    imgType=4;
                    numROI=size(seeds.vol,4);
                end
                fprintf('%d unique ROI indices identified in %dD volume \n%s\n',numROI,imgType,fileMasked)

                % read in label file
                if exist(paths.DWI.seedLabels,'file')
                    fid=fopen(paths.DWI.seedLabels);
                    slabels=textscan(fid,'%s');
                    fclose(fid);
                    if length(slabels{1})~=numROI
                        disp('Number of labels does not match number of ROI.')
                        disp('Double check the label text file for proper size/format.')
                        disp('......')
                        disp('Proceeding without labels. ROI will be identified by integer values')
                        slabelON=0;
                    else
                        slabelON=1;
                    end
                else
                    disp('Seed label file does not exist. Proceeding without labels.')
                    slabelON=0;
                end

                for s=1:numROI
                    if slabelON==1
                        fileOut=fullfile(paths.DWI.nodeTracks,sprintf('seed_MNI_%s.nii.gz',slabels{1}{s}));
                    elseif slabelON==0
                        if imgType==3
                            fileOut=fullfile(paths.DWI.nodeTracks,sprintf('seed_MNI_%d.nii.gz',seedIds(s)));
                        elseif imgType==4
                            fileOut=fullfile(paths.DWI.nodeTracks,sprintf('seed_MNI_%d.nii.gz',s));
                        end
                    end
                    %isolate seed ROI image
                    if imgType==3
                        sentence=sprintf('%s/fslmaths %s -thr %d -uthr %d %s',...
                            paths.FSL,fileMasked,seedIds(s),seedIds(s),fileOut);
                        [~,result]=system(sentence);
                    elseif imgType==4
                        sentence=sprintf('%s/fslroi %s %s %d 1',...
                            paths.FSL,fileMasked,fileOut,numROI);
                        [~,result]=system(sentence);
                    end
                    seedFiles{s,1}=fileOut;
                end
            else
                warning('MNI Seed image not found. Exiting ...')
                return
            end

        otherwise
            warning('Not a recognized input type for params.DWI.seedType. Exiting...')
            return
    end

    % Parse target inputs
    if length(params.DWI.SeedType)==2
        switch params.DWI.SeedType(2)
            case 1 % node targets
            disp('-------------------------------')
            disp('   Working on: Nodes as Targets')
            disp('-------------------------------')
            if ~isempty(configs.DWI.TargetIdx) 
                if exist(paths.DWI.targetLabels,'file')
                    % read in label file
                    fid=fopen(paths.DWI.targetLabels);
                    tlabels=textscan(fid,'%s');
                    fclose(fid);
                    if length(tlabels{1})~=length(configs.DWI.TargetIdx)
                        disp('Number of labels does not match number of ROI.')
                        disp('Double check the label text file for proper size/format.')
                        disp('......')
                        disp('Proceeding without labels. ROI will be identified by integer values')
                        tlabelON=0;
                    else
                        tlabelON=1;
                    end
                else
                    disp('Target label file does not exist. Proceeding without labels.')
                    tlabelON=0;
                end
                for s=1:length(configs.DWI.TargetIdx) % for every seed in directory
                    % generate seed volume
                    for k=1:length(parcs.pnodal)
                        if parcs.pnodal(k).true == 1
                            fileIn=fullfile(paths.DWI.dir,sprintf('rT1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
                            if tlabelON==1
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('target_%s_%s.nii.gz',parcs.plabel(k).name,tlabels{1}{s}));
                            elseif tlabelON==0
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('target_%s_%d.nii.gz',parcs.plabel(k).name,configs.DWI.TargetIdx(s)));
                            end
                            sentence=sprintf('%s/fslmaths %s -thr %d -uthr %d %s',...
                                paths.FSL,fileIn,configs.DWI.TargetIdx(s),configs.DWI.TargetIdx(s),fileOut);
                            [~,result]=system(sentence);
                            % Store names of seed files for probtrackX use
                            targetFiles{s,1}=fileOut;
                        end
                    end
                end
            else
                warning('configs.DWI.TargetIdx is empty. Exiting...')
                return
            end  

            case 2 % MNI targets
                disp('-------------------------------')
                disp('        Working on: MNI Targets')
                disp('-------------------------------')
                if exist(paths.DWI.MNItargetImage,'file')
                    % Prepare the target images
                    % Tranform to diffusion space
                    %concatenate linear transformations form T1
                    MNI2T1_linear=fullfile(paths.T1.reg,'MNI2T1_linear.mat');
                    if ~exist(MNI2T1_linear,'file')
                        dof6_2T1=fullfile(paths.T1.reg,'MNI2T1_dof6.mat');
                        dof12_2T1=fullfile(paths.T1.reg,'MNI2T1_dof12.mat');
                        sentence=sprintf('%s/convert_xfm -omat %s -concat %s %s',...
                            paths.FSL,MNI2T1_linear,dof12_2T1,dof6_2T1);
                        [~,result]=system(sentence);
                    end
                    % apply MNI2T1 transformation
                    MNI2T1_warp=fullfile(paths.T1.reg,'MNI2T1_warp.nii.gz');
                    fileT1ref=fullfile(paths.T1.dir,'T1_fov_denoised.nii');
                    [~,name,ext]=fileparts(paths.DWI.MNItargetImage);
                    T1targetImage=fullfile(paths.DWI.nodeTracks,sprintf('T1_%s',[name ext]));
                    sentence=sprintf('%s/applywarp -i %s -r %s -o %s -w %s --postmat=%s --interp=nn',...
                        paths.FSL,paths.DWI.MNItargetImage,fileT1ref,T1targetImage,MNI2T1_warp,MNI2T1_linear);
                    [~,result]=system(sentence);
                    % apply T1 to DWI transformation
                    fileMat2FA=fullfile(paths.DWI.dir,'T1_2FAdof6.mat');
                    filewarp=fullfile(paths.DWI.dir,'T12FA_warpfield.nii.gz');
                    DWItargetImage=fullfile(paths.DWI.nodeTracks,sprintf('DWI_%s',[name ext]));
                    fileRef=fullfile(paths.DWI.DTIfit,'3_DWI_FA.nii.gz');
                    if exist(filewarp,'file')
                        sentence=sprintf('%s/applywarp -i %s -o %s -r %s -w %s --interp=nn',...
                            paths.FSL,T1targetImage,DWItargetImage,fileRef,filewarp);
                        [~,result]=system(sentence);
                    else
                        sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                            paths.FSL,T1targetImage,fileRef,fileMat2FA,DWItargetImage);
                        [~,result]=system(sentence);
                    end

                    % mask & dilate image within gray matter
                    GM_Mask=fullfile(paths.DWI.dir,'rT1_GM_mask.nii.gz');
                    fileMasked=fullfile(paths.DWI.nodeTracks,sprintf('DWI_GM_%s',[name ext]));
                    sentence=sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,DWItargetImage,GM_Mask,fileMasked);
                    [~,result]=system(sentence);
%                     sentence=sprintf('%s/fslmaths %s -dilD -mas %s %s',paths.FSL,fileMasked,GM_Mask,fileMasked);
%                     [~,result]=system(sentence);

                    % Read in and sort target image
                    targets=MRIread(fileMasked);
                    num4rthDIM=size(targets.vol,4);
                    if num4rthDIM==1 %3D
                        imgType=3;
                        targetIds=unique(targets.vol);
                        targetIds(targetIds==0)=[];
                        numROI=length(targetIds);   
                    elseif num4rthDIM > 1 %4D
                        imgType=4;
                        numROI=size(targets.vol,4);
                    end
                    fprintf('%d unique ROI indices identified in %dD volume \n%s\n',numROI,imgType,fileMasked)

                    % read in label file
                    if exist(paths.DWI.targetLabels,'file')
                        fid=fopen(paths.DWI.targetLabels);
                        tlabels=textscan(fid,'%s');
                        fclose(fid);
                        if length(tlabels{1})~=numROI
                            disp('Number of labels does not match number of ROI.')
                            disp('Double check the label text file for proper size/format.')
                            disp('......')
                            disp('Proceeding without labels. ROI will be identified by integer values')
                            tlabelON=0;
                        else
                            tlabelON=1;
                        end
                    else
                        disp('Target label file does not exist. Proceeding without labels.')
                        tlabelON=0;
                    end

                    for s=1:numROI
                        if tlabelON==1
                            fileOut=fullfile(paths.DWI.nodeTracks,sprintf('target_MNI_%s.nii.gz',tlabels{1}{s}));
                        elseif tlabelON==0
                            if imgType==3
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('target_MNI_%d.nii.gz',targetIds(s)));
                            elseif imgType==4
                                fileOut=fullfile(paths.DWI.nodeTracks,sprintf('target_MNI_%d.nii.gz',s));
                            end
                        end
                        %isolate target ROI image
                        if imgType==3
                            sentence=sprintf('%s/fslmaths %s -thr %d -uthr %d %s',...
                                paths.FSL,fileMasked,targetIds(s),targetIds(s),fileOut);
                            [~,result]=system(sentence);
                        elseif imgType==4
                            sentence=sprintf('%s/fslroi %s %s %d 1',...
                                paths.FSL,fileMasked,fileOut,numROI);
                            [~,result]=system(sentence);
                        end
                        targetFiles{s,1}=fileOut;
                    end
                else
                    warning('MNI Target image not found. Exiting ...')
                    return
                end
            otherwise
                warning('Not a recognized input type for params.DWI.seedType. Exiting...')
                return
        end
    end

    % set up inputs in pieces
    sampleBaseName=fullfile(paths.DWI.bpxOUT,'merged');
    brainmaskVol=fullfile(paths.DWI.bpxOUT,'nodif_brain_mask.nii.gz');
    hardset='-l --modeuler --onewaycondition --fibthresh=0.01 --forcedir --opd --ompl';
    userset=sprintf('-c %.2f -S %d --steplength=%.2f -P %d --distthresh=%.2f --sampvox=%.2f',...
        configs.DWI.curvature,configs.DWI.maxSteps,configs.DWI.stepLength,configs.DWI.numSamples,configs.DWI.disttresh,configs.DWI.sampvox);

    if params.DWI.network == 1
        outDir=fullfile(paths.DWI.nodeTracks,configs.DWI.netName);
        if ~exist(outDir,'dir')
            mkdir(outDir)
        end
        fileMasks=fullfile(outDir,'masks.txt');
        fileID=fopen(fileMasks,'w');
        for m=1:length(seedFiles)
            fprintf(fileID,'%s\n',seedFiles{m,1});
        end
        for m=1:length(targetFiles)
            fprintf(fileID,'%s\n',targetFiles{m,1});
        end
        fclose(fileID);
        codeset=sprintf('--network -x %s -s %s -m %s --dir=%s',fileMasks,sampleBaseName,brainmaskVol,outDir);
    else
        for i=1:length(seedFiles)
            [~,Sname,~]=fileparts(seedFiles{i,1}); %remove path and .gz
            [~,Sname,~]=fileparts(Sname); % removes .nii
            if exist('targetFiles','var')
                for j=1:length(targetFiles)
                    [~,Tname,~]=fileparts(targetFiles{j,1}); %remove path and .gz
                    [~,Tname,~]=fileparts(Tname); % removes .nii
                    outDir=fullfile(paths.DWI.nodeTracks,sprintf('%s-%s',Sname,Tname));
                    if ~exist(outDir,'dir')
                     mkdir(outDir)
                    end
                    %write seed/target to mask file
                    fileMasks=fullfile(outDir,'targets.txt');
                    fileID=fopen(fileMasks,'w');
                    fprintf(fileID,'%s\n',targetFiles{j,1});
                    fclose(fileID);
                    codeset=sprintf('-x %s -s %s -m %s --dir=%s --targetmasks=%s --os2t --otargetpaths',seedFiles{i,1},sampleBaseName,brainmaskVol,outDir,fileMasks);
                    % piece together input & run
                    sentence=sprintf('%s/probtrackx2 %s %s %s',paths.FSL,codeset,hardset,userset);
                    [~,result]=system(sentence);
                    disp(result)
                end
            else
                codeset=sprintf('-x %s.nii.gz -s %s -m %s --dir=%s',seedFiles{i,1},sampleBaseName,brainmaskVol,outDir);
                % piece together input & run
                sentence=sprintf('%s/probtrackx2 %s %s %s',paths.FSL,codeset,hardset,userset);
                [~,result]=system(sentence);
                disp(result) 
            end
        end
    end
    % piece together input & run
    sentence=sprintf('%s/probtrackx2 %s %s %s',paths.FSL,codeset,hardset,userset);
    [~,result]=system(sentence);
    disp(result)  
end
