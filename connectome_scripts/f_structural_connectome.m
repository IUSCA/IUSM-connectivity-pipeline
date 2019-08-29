function [paths,flags,configs,parcs]=f_structural_connectome(paths,flags,configs,parcs)
%                       F_STRUCTURAL_CONNECTOME
% Does the registration of anatomical and parcellation data, refrorms
% tractography, and generates connectivity matrices of sctructural data.
%
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

%% registration of B0 to T1
if flags.DWI.reg2T1==1
    disp('------------------------')
    disp('Registration of B0 to T1')
    disp('------------------------')
    
    % rigid body of bo to T1
    disp('rigid body dof 6 to T1')
    fileIn =  fullfile(paths.DWI.DTIfit,'b0_1st.nii.gz');
    fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
    fileMat1 = fullfile(paths.DWI.dir,'DWI_2_T1_dof6.mat');
    fileOut = fullfile(paths.DWI.dir,'b0_avg.T1.dof6.nii.gz');
    sentence = sprintf('%s/flirt -in %s -ref %s -omat %s -cost normmi -dof 6 -interp spline -out %s',...
        paths.FSL,fileIn,fileRef,fileMat1,fileOut);
    [~,result] = system(sentence); %#ok<*ASGLU>
    
    % bbr rigid body of bo to T1
    disp('rigid body dof bbr to T1')
    fileIn = fullfile(paths.DWI.dir,'b0_avg.T1.dof6.nii.gz');
    fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
    fileWMmask = fullfile(paths.T1.dir,'T1_WM_mask.nii.gz');
    fileMat2 = fullfile(paths.DWI.dir,'DWI_2_T1_bbr.mat');
    sentence = sprintf('%s/flirt -in %s -ref %s -wmseg %s -cost bbr -omat %s',...
        paths.FSL,fileIn,fileRef,fileWMmask,fileMat2);
    [~,result] = system(sentence);

    fileMat1 = fullfile(paths.DWI.dir,'DWI_2_T1_dof6.mat');
    fileMat2 = fullfile(paths.DWI.dir,'DWI_2_T1_bbr.mat');
    disp('concatenating matrices')
    fileMatJoint = fullfile(paths.DWI.dir,'DWI_2_T1_joint.mat');
    sentence = sprintf('%s/convert_xfm -omat %s -concat %s %s',...
        paths.FSL,fileMatJoint,fileMat2,fileMat1);
    [~,result] = system(sentence);
    
    % udpate bvec and copy bval
    fileMat = fullfile(paths.DWI.dir,'DWI_2_T1_joint.mat');    
    matrix = dlmread(fileMat);    
    bvecs =dlmread(fullfile(paths.DWI.EDDY,'eddy_output.eddy_rotated_bvecs'));
    bvecs_corr = f_correct_bvec_from_matrix(bvecs,matrix);
 
    % apply DWI_2_T1_joint.mat when possible
    if any((sum(bvecs_corr(2:end,:).^2,2)-1)>0.05)
        fprintf('bvecs not unitary, apply only dof6 transformation \n');
        fileMat = fullfile(paths.DWI.dir,'DWI_2_T1_dof6.mat');
        matrix = dlmread(fileMat);    
        bvecs =dlmread(fullfile(paths.DWI.EDDY,'eddy_output.eddy_rotated_bvecs'));
        bvecs_corr = f_correct_bvec_from_matrix(bvecs,matrix);
        if any((sum(bvecs_corr(2:end,:).^2,2)-1)>0.05)
            fprintf('bvecs not unitary even with dof6 only..exiting \n');
            return;
        else
            dlmwrite(fullfile(paths.DWI.dir,'4_DWI.bvec'),bvecs_corr,'delimiter',' ','precision','%.6f');
            
            disp('applying transformation')
            fileIn = fullfile(paths.DWI.EDDY,'eddy_output.nii.gz');
            fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
            fileOut = fullfile(paths.DWI.dir,'4_DWI.nii.gz');
            sentence = sprintf('%s/flirt -nosearch -in %s -ref %s -applyxfm -init %s -out %s -interp spline',...
                paths.FSL,fileIn,fileRef,fileMat,fileOut);
            [~,result] = system(sentence);
        end
    else
    dlmwrite(fullfile(paths.DWI.dir,'4_DWI.bvec'),bvecs_corr,'delimiter',' ','precision','%.6f');
    
        disp('applying transformation')
        fileIn = fullfile(paths.DWI.EDDY,'eddy_output.nii.gz');
        fileRef = fullfile(paths.T1.dir,'T1_fov_denoised.nii');
        fileOut = fullfile(paths.DWI.dir,'4_DWI.nii.gz');
        sentence = sprintf('%s/flirt -nosearch -in %s -ref %s -applyxfm -init %s -out %s -interp spline',...
                paths.FSL,fileIn,fileRef,fileMat,fileOut);
        [~,result] = system(sentence);
    end
    
    fileIn = fullfile(paths.DWI.dir,'4_DWI.nii.gz');
    fileOut = fullfile(paths.DWI.dir,'4_avg_brain.nii.gz');
    sentence = sprintf('%s/fslroi %s %s 0 1',paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.DWI.dir,'4_avg_brain.nii.gz');
    sentence = sprintf('%s/bet %s %s -f 0.20 -g 0 -m -t',paths.FSL,fileIn,fileIn);
    [~,result] = system(sentence);
    fileIn = fullfile(paths.DWI.dir,'4_avg_brain_mask.nii.gz');
    sentence =sprintf('%s/fslmaths %s -fillh %s',paths.FSL,fileIn,fileIn);
    [~,result]=system(sentence);
end

%% Generation Tissue Masks (seed; fiber; GM_WM boundary)
if flags.DWI.tissueMasks==1
    disp('------------')
    disp('Tissue Masks')
    disp('------------')
    
   % WM ring
    fileIn = fullfile(paths.T1.dir,'T1_GM_mask');
    fileOut = fullfile(paths.T1.dir,'T1_GM_mask_dil');
    sentence = sprintf('%s/fslmaths %s -dilD %s',...
        paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);

    fileIn = fullfile(paths.T1.dir,'T1_GM_mask');
    fileAdd = fullfile(paths.T1.dir,'T1_WM_mask');
    fileOut = fullfile(paths.T1.dir,'T1_GMWM_mask');
    sentence = sprintf('%s/fslmaths %s -add %s %s',paths.FSL,fileIn,fileAdd,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.T1.dir,'T1_GM_mask_dil');
    fileMul = fullfile(paths.T1.dir,'T1_WM_mask');
    fileOut = fullfile(paths.T1.dir,'WM_ring');
    sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn,fileMul,fileOut);
    [~,result] = system(sentence);
    
    % GM ring
    fileIn = fullfile(paths.T1.dir,'T1_WM_mask');
    fileOut = fullfile(paths.T1.dir,'T1_WM_dil');
    sentence = sprintf('%s/fslmaths %s -dilD %s',...
        paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.T1.dir,'T1_WM_dil');
    fileMul = fullfile(paths.T1.dir,'T1_GM_mask');
    fileOut = fullfile(paths.T1.dir,'GM_ring');
    sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn,fileMul,fileOut);
    [~,result] = system(sentence);
    
    % mask 4 fibers
    fileIn = fullfile(paths.T1.dir,'T1_WM_mask');
    fileAdd = fullfile(paths.T1.dir,'GM_ring');
    fileOut = fullfile(paths.DWI.dir,'mask4fibers');
    sentence = sprintf('%s/fslmaths %s -add %s %s',paths.FSL,fileIn,fileAdd,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.DWI.dir,'mask4fibers');
    fileOut = fullfile(paths.DWI.dir,'mask4fibers');
    sentence = sprintf('%s/fslmaths %s -bin %s',paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.DWI.dir,'mask4fibers');
    fileMul = fullfile(paths.DWI.dir,'4_avg_brain_mask');
    fileOut = fullfile(paths.DWI.dir,'mask4fibers');
    sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn,fileMul,fileOut);
    [~,result] = system(sentence);
    
    % mask 4 seeds
    fileIn1 = fullfile(paths.T1.dir,'GM_ring');
    fileIn2 = fullfile(paths.T1.dir,'WM_ring');
    fileOut = fullfile(paths.DWI.dir,'mask4seeds');
    sentence = sprintf('%s/fslmaths %s -add %s %s',paths.FSL,fileIn1,fileIn2,fileOut);
    [~,result] = system(sentence);
    sentence = sprintf('%s/fslmaths %s -bin %s',paths.FSL,fileOut,fileOut);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.DWI.dir,'mask4seeds');
    fileMul = fullfile(paths.DWI.dir,'4_avg_brain_mask');
    fileOut = fullfile(paths.DWI.dir,'mask4seeds');
    sentence = sprintf('%s/fslmaths %s -mul %s %s',paths.FSL,fileIn,fileMul,fileOut);
    [~,result] = system(sentence);
end

%% Run Camino tensor modeling and tractography.
if flags.DWI.Camino==1
caminoMem = sprintf('export CAMINO_HEAP_SIZE=%d',configs.DWI.HeapSizeCamino);
   
    paths.DWI.Camino = fullfile(paths.DWI.dir,'Camino');
    if configs.DWI.CaminoReset==1
        if(exist(paths.DWI.Camino,'dir'))
            sentence = sprintf('rm -rf %s',paths.DWI.Camino);
            [~,result] = system(sentence);        
        end
        mkdir(paths.DWI.Camino);
    end
    if(~exist(paths.DWI.Camino,'dir'))
        mkdir(paths.DWI.Camino);
    end

%% camino processing
    if flags.DWI.getTensor==1
    disp('-----------------')
    disp('Camino Processing')
    disp('-----------------')
    
    % generate a scheme file for camino to use with direction orientaitons
    fileDirs = fullfile(paths.DWI.dir,'4_DWI.bvec');
    fileBvals = fullfile(paths.DWI.dir,'0_DWI.bval');
    fileScheme = fullfile(paths.DWI.Camino,'schemeFileCamino.txt');
    sentence = sprintf('fsl2scheme -bvecfile %s -bvalfile %s -flipz -bscale 1 > %s',...
        fileDirs,fileBvals,fileScheme);
    [~,result] = system(sentence);
    disp(result)
    DTIfile = fullfile(paths.DWI.dir,'4_DWI.nii.gz');

    if ~exist(DTIfile,'file')
        error('DTIfile not found for subj!')
    end
    if ~exist(fileScheme,'file')
        error('Scheme file for Camino not found for subj!')
    end
    maskBrain = fullfile(paths.T1.dir,'T1_GMWM_mask.nii.gz');
    mask4Fibers = fullfile(paths.DWI.dir,'mask4fibers.nii.gz');

    % convert DWI from image to voxel data.
    dwiFile = fullfile(paths.DWI.Camino,'dwi.Bfloat');
    sentence = sprintf('%s;image2voxel -4dimage %s -outputfile %s',caminoMem,DTIfile,dwiFile);
    [~,result] = system(sentence);
    % Fit a single tensor model
        % model : nonlinear optimization, constrained to be positive semi-definite.
    dtModel = fullfile(paths.DWI.Camino,'dt_nonlinear.Bdouble');
    sentence = sprintf('%s;modelfit -model nldt_pos -schemefile %s -inputfile %s -bgmask %s > %s',caminoMem,fileScheme,dwiFile,fullfile(paths.DWI.dir,'4_avg_brain_mask.nii.gz'),dtModel);
    [~,result] = system(sentence);
    disp(result)
    % Extract a fractional anisotropy and mean diffusivity image.
    DTImaps = ['fa';'md'];
    for i=1:length(DTImaps(:,1))
    Qntfile = fullfile(paths.DWI.Camino,DTImaps(i,:));
    sentence = sprintf('%s;cat %s | %s | voxel2image -outputroot %s -header %s',caminoMem,dtModel,DTImaps(i,:),Qntfile,maskBrain);
    [~,result] = system(sentence);
    disp(result)
    end
    % Extract other diffusivity maps
        % Convert data to nifti
        TensorOut = fullfile(paths.DWI.Camino,'tensor_');
        sentence = sprintf('%s;dt2nii -inputfile %s -inputdatatype double -header %s -outputroot %s',caminoMem,dtModel,maskBrain,TensorOut);
        [~,result] = system(sentence);
        
        % Extract radial diffusivity
        TensorImg = fullfile(paths.DWI.Camino,'tensor_dt.nii.gz');
        RDout = fullfile(paths.DWI.Camino,'rd.nii.gz');
        sentence = sprintf('TVtool -in %s -rd -out %s',TensorImg,RDout);
        [~,result] = system(sentence);
        % Extract axial diffusivity
        ADout = fullfile(paths.DWI.Camino,'ad.nii.gz');
        sentence = sprintf('TVtool -in %s -ad -out %s',TensorImg,ADout);
        [~,result] = system(sentence);

    % Classify voxels as isotropic, single tensor, or crossing fibre
    fileVC = fullfile(paths.DWI.Camino,sprintf('dwiVC_%d_%d.Bfloat',configs.DWI.order2Threshold,configs.DWI.order4Threshold));

    sentence = sprintf('%s;voxelclassify -schemefile %s -bgmask %s -order 4 -ftest 1.0E-%d 1.0E-%d 1.0E-%d < %s > %s',...
        caminoMem,fileScheme,mask4Fibers,configs.DWI.order2Threshold,configs.DWI.order4Threshold,configs.DWI.order6Threshold,dwiFile,fileVC);
    display(sentence);
    [~,result] = system(sentence);
    disp(result)
    
    fileVCnii = fileVC(1:end-7); %remove suffix ".Bfloat" from the name

    sentence = sprintf('cat %s | voxel2image -inputdatatype int -header %s -outputroot %s',fileVC,maskBrain,fileVCnii);
    display(sentence)
    [~,result] = system(sentence);
    disp(result)

    vol_VC = MRIread(sprintf('%s.nii.gz',fileVCnii));
    VCcomps = bwconncomp(vol_VC.vol>=4); % get compononets for voxels classified as multi-tensor.
    for icomps = 1:VCcomps.NumObjects
        comp = VCcomps.PixelIdxList{icomps};
        if length(comp)<configs.DWI.clusterMinSizeMultiTensor % if the cluster of multi-tensor voxels do not have more than 16
            vol_VC.vol(comp) = 2; %degrade from multi-tensor to mono-tensor
        end                    
    end
    fileVCnii = sprintf('%s_clean.nii.gz',fileVCnii);
    MRIwrite(vol_VC,fileVCnii);

    fileVCclean = fullfile(paths.DWI.Camino,sprintf('dwiVC_%d_%d_clean.Bint',configs.DWI.order2Threshold,configs.DWI.order4Threshold));
    fileImageList = fullfile(paths.DWI.Camino,'vcvol.txt');
    fid = fopen(fileImageList,'wt');
    fprintf(fid,'%s\n',sprintf('dwiVC_%d_%d_clean.nii.gz',configs.DWI.order2Threshold,configs.DWI.order4Threshold));
    fclose(fid);
    sentence = sprintf('%s;image2voxel -imagelist %s -outputdatatype int > %s',caminoMem,fileImageList,fileVCclean);
    display(sentence);
    [~,result] = system(sentence);

    multitensorModel = fullfile(paths.DWI.Camino,'multitensor.Bdouble');
    sentence = sprintf('%s;multitenfit -inputfile %s -bgmask %s -schemefile %s -voxclassmap %s > %s',...
        caminoMem,dwiFile,mask4Fibers,fileScheme,fileVCclean,multitensorModel);
    display(sentence);
    [~,result] = system(sentence);

    eigFile = fullfile(paths.DWI.Camino,'multitensor_eigs.Bdouble');
    sentence = sprintf('%s;dteig -inputmodel multitensor -maxcomponents 2 < %s > %s',...
    caminoMem,multitensorModel,eigFile);
    disp(sentence); 
    [~,result] = system(sentence);
    end
    
%% deterministic multitensor tractography (2-tensor max)
    if flags.DWI.Deterministic==1
    disp('--------------------------')
    disp('Deterministic Tractography')
    disp('--------------------------')   
    
    multitensorModel = fullfile(paths.DWI.Camino,'multitensor.Bdouble');
    if ~exist(multitensorModel,'file')
        warning('%s not exist; Check that getTensor flag was ran.',multitensorModel)
        return
    else
    
    paths.DWI.Dtrack = fullfile(paths.DWI.Camino,'Determ_Tractography');
    
    if(exist(paths.DWI.Dtrack,'dir'))
        sentence = sprintf('rm -rf %s',paths.DWI.Dtrack);
        [~,result] = system(sentence);        
    end
        mkdir(paths.DWI.Dtrack);
    if(~exist(paths.DWI.Dtrack,'dir'))
        mkdir(paths.DWI.Dtrack);
    end
    
    mask4fibers = fullfile(paths.DWI.dir,'mask4fibers.nii.gz');
    if configs.DWI.seedsWMborder==1 %seeds only in the interface of WM and GM (THE RING)
        mask4seeds = fullfile(paths.DWI.dir,'mask4seeds.nii.gz');
    else
        mask4seeds = fullfile(paths.T1.dir,'WM_ring.nii.gz');
    end
    
    fileVCclean = fullfile(paths.DWI.Camino,sprintf('dwiVC_%d_%d_clean.Bint',...
        configs.DWI.order2Threshold,configs.DWI.order4Threshold));
    fileTracts = fullfile(paths.DWI.Dtrack,'allTractsMultiTensor.Bfloat');
    sentence = sprintf('%s;track -inputmodel multitensor -maxcomponents 2 -voxclassmap %s -inputfile %s -seedfile %s -curvethresh %d -curveinterval 5 -anisfile %s -anisthresh %s -stepsize %0.2f -iterations 1 -outputfile %s',...
        caminoMem,fileVCclean,multitensorModel,mask4seeds,configs.DWI.CURVthresh,mask4fibers,0.5,configs.DWI.stepSize,fileTracts);
    disp(sentence);
    [~,~] = system(sentence);
    end
    end

%% Create trk streamline files.
    if flags.DWI.FiberFiles==1
    disp('------------------')
    disp('Create Fiber Files')
    disp('------------------')

    paths.DWI.Dtrack = fullfile(paths.DWI.Camino,'Determ_Tractography');

    if exist(fullfile(paths.DWI.Camino,'Determ_Tractography'),'dir')
        fileTracts =  fullfile(paths.DWI.Dtrack,'allTractsMultiTensor.Bfloat');
        fileTractsFltrd = fullfile(paths.DWI.Dtrack,'allTractsMultiTensorFltrd.Bfloat');
        fileHeader = fullfile(paths.DWI.dir,'4_avg_brain.nii.gz');
    end

    % filter fiber length    
    sentence = sprintf('%s;procstreamlines -mintractlength %d -maxtractlength %d -header %s <%s >%s',...
        caminoMem,configs.DWI.LengthMin,configs.DWI.LengthMax,fileHeader,fileTracts,fileTractsFltrd);
    disp(sentence)
    [~,~] = system(sentence);

    % compute fiber-length
    mask4fibers = fullfile(paths.DWI.dir,'mask4fibers.nii.gz');
    fileTractLength = fullfile(paths.DWI.Dtrack,'tractLength.Bdouble');
    sentence = sprintf('cat %s | tractstats -header %s -tractstat length > %s',fileTractsFltrd,mask4fibers,fileTractLength);
    disp(sentence);
    [~,result] = system(sentence);
    disp(result)
    
    mask4Fibers = fullfile(paths.DWI.dir,'mask4fibers.nii.gz');
    FAfile = fullfile(paths.DWI.Camino,'fa.nii.gz');
    
    volFA = MRIread(FAfile);
    volMASK4fibers = MRIread(mask4Fibers);
    volFA4fibers = volFA;
    volFA4fibers.vol = volFA.vol.*volMASK4fibers.vol.*(volFA.vol>configs.DWI.FAthresh);
    FAfile4fibers = fullfile(paths.DWI.Camino,'fa4fibers.nii.gz');
    MRIwrite(volFA4fibers,FAfile4fibers,'double');
    clear volFA volFA4fibers;

    % convertion from vtk (Camino format) to trk (trackVis format)
    % set tmp files
    params.subjName=extractAfter(paths.subject,sprintf('%s/',paths.data));
    tmp_tracts = sprintf('/tmp/%s_tracts.Bfloat',params.subjName);
    tmp_trk = sprintf('/tmp/%s_tracts.trk',params.subjName);
    tmp_vtk = sprintf('/tmp/%s_tracts.vtk',params.subjName);
    tmp_fa = sprintf('/tmp/%s_fa.nii.gz',params.subjName);
    
    % copy tracts file to tmp
    sentence = sprintf('cp %s %s',fileTractsFltrd,tmp_tracts);
    disp(sentence)
    [~,result] = system(sentence);
    
    % copy fa4fibers file to tmp
    fileFA = fullfile(paths.DWI.Camino,'fa4fibers.nii.gz');
    sentence = sprintf('cp %s %s',fileFA,tmp_fa);
    disp(sentence)
    [~,result] = system(sentence);
    
    sentence = sprintf('%s;camino_to_trackvis -i %s -o %s --nifti %s --phys-coords',paths.CamTrackSetup,tmp_tracts,tmp_trk,tmp_fa);
    disp(sentence)
    [~,result] = system(sentence);
    disp(result)
    
    % copy trk from tmp to subject-folder
    fileTRKphys = fullfile(paths.DWI.Dtrack,'3D_fibers_phys.trk');
    sentence = sprintf('cp %s %s',tmp_trk,fileTRKphys);
    disp(sentence)
    [~,result] = system(sentence);

    % convertion to vtk
    sentence = sprintf('%s;vtkstreamlines < %s > %s',caminoMem,tmp_tracts,tmp_vtk);
    disp(sentence)
    [~,result] = system(sentence);
    
    % copy vtk from tmp to subject-folder
    fileVTK = fullfile(paths.DWI.Dtrack,'3D_fibers.vtk');
    sentence = sprintf('cp %s %s',tmp_vtk,fileVTK);
    disp(sentence)
    [~,result] = system(sentence);
    
    % clean tmp files
    sentence = sprintf('rm -f %s %s %s %s',tmp_tracts,tmp_trk,tmp_vtk,tmp_fa);
    disp(sentence)
    [~,result] = system(sentence);  
    end
  
%% Generate tract-based metric matrices.
    if flags.DWI.genMats==1
    disp('-----------------')
    disp('Generate Matrices')
    disp('-----------------')
    
    paths.DWI.Dtrack = fullfile(paths.DWI.Camino,'Determ_Tractography');
    if exist(fullfile(paths.DWI.Camino,'Determ_Tractography'),'dir')
      fileTractsFltrd = fullfile(paths.DWI.Dtrack,'allTractsMultiTensorFltrd.Bfloat');
    end
    paths.DWI.matrices = fullfile(paths.DWI.Dtrack,'matrices');
    if ~exist(paths.DWI.matrices,'dir')
        mkdir(paths.DWI.matrices);
    end

    FAfile = fullfile(paths.DWI.Camino,'fa.nii.gz');
    
    for k=1:length(parcs.pdir)
        if parcs.pnodal(k).true == 1    
            suffix = sprintf('%s_multiTensor',parcs.plabel(k).name);
            maskGMparc = fullfile(paths.T1.dir,sprintf('T1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
            MATRIXfile = fullfile(paths.DWI.matrices,sprintf('M_%s_',suffix));
    
    %% number of fibers and <FA>
            sentence = sprintf('%s; conmat -inputfile %s -targetfile %s -scalarfile %s -tractstat mean -outputroot %s'...
                ,paths.CaminoSetup,fileTractsFltrd,maskGMparc,FAfile,MATRIXfile);
            
            [~,result] = system(sentence);
            fileMnf = fullfile(paths.DWI.matrices,sprintf('M_%s_sc.csv',suffix));
            Mnf = csvread(fileMnf,1,0);
            fileMfaAvg = fullfile(paths.DWI.matrices,sprintf('M_%s_ts.csv',suffix));
            MfaAvg = csvread(fileMfaAvg,1,0); %#ok<*NASGU>
    %% minimum FA
            sentence = sprintf('%s; conmat -inputfile %s -targetfile %s -scalarfile %s -tractstat min -outputroot %s'...
                ,paths.CaminoSetup,fileTractsFltrd,maskGMparc,FAfile,MATRIXfile);
            [~,result]=system(sentence);
            fileMfaMin = fullfile(paths.DWI.matrices,sprintf('M_%s_ts.csv',suffix));
            MfaMin = csvread(fileMfaMin,1,0);

    %% longest length
    sentence = sprintf('%s; conmat -inputfile %s -targetfile %s -tractstat length -outputroot %s'...
        ,paths.CaminoSetup,fileTractsFltrd,maskGMparc,MATRIXfile);
            [~,result]=system(sentence);
            fileMll = fullfile(paths.DWI.matrices,sprintf('M_%s_ts.csv',suffix));
            Mll = csvread(fileMll,1,0);

    %% weighted density and mean surface area
            GMring_file = fullfile(paths.T1.dir,'GM_ring.nii.gz');
            GMring_vol = MRIread(GMring_file); GMring_vol = GMring_vol.vol;
            GMparc_vol = MRIread(maskGMparc);
            GMparc_vol = GMparc_vol.vol;
            GMparc_ring_vol = GMparc_vol .* GMring_vol;
    
            numROIs = size(Mnf,1);
            Sroi = zeros(numROIs,1);
                for roi_index=1:numROIs
                    Sroi(roi_index) = nnz(GMparc_ring_vol==roi_index);
                end
            Msurf = zeros(numROIs,numROIs);
                for roi_i=1:numROIs
                    for roi_j=1:numROIs
                       Msurf(roi_i,roi_j) = (Sroi(roi_i) + Sroi(roi_j))/2;
                    end
                end
            Mw = Mnf./Msurf;

%% Save matrices
            fname = fullfile(paths.DWI.matrices,sprintf('Mnf_%s.mat',suffix));
            varname = 'Mnf';
            save(fname,varname); 

            fname =  fullfile(paths.DWI.matrices,sprintf('MfaAvg_%s.mat',suffix));
            varname = 'MfaAvg';
            save(fname,varname); 

            fname =  fullfile(paths.DWI.matrices,sprintf('MfaMin_%s.mat',suffix));
            varname = 'MfaMin';
            save(fname,varname); 

            fname =  fullfile(paths.DWI.matrices,sprintf('Mll_%s.mat',suffix));
            varname = 'Mll';
            save(fname,varname);

            fname =  fullfile(paths.DWI.matrices,sprintf('Msurf_%s.mat',suffix));
            varname = 'Msurf';
            save(fname,varname);

            fname =  fullfile(paths.DWI.matrices,sprintf('Mw_%s.mat',suffix));
            varname = 'Mw';
            save(fname,varname);
        end
    end
    end
end