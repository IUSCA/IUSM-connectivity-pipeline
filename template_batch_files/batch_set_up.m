%                               BATCH_SET_UP
% This contains all flags and configs required by the pipeline to process
% the data. You may edit this as necessary depending on what portions of
%                       the pipeline you wish to run.
%
%      Evgeny Chumin, Indiana University School of Medicine, 2018
%%
                  %-------------------------------%
                  %  SET UP PARCELLATION OPTIONS  %
                  %-------------------------------%
% IMPORTANT: If you want to introduce a new parcellation into the pipeline,
%   follow these steps and refer to existing parcellations as examples.
%   1. In connectome_scripts/templates/MNIparcs:
%       -Create a directory with the same name as the parcellation (just
%       without the .nii.gz)
%       -In that directory place the .nii.gz parcellation image that is in
%       MNI152 space (e.g. registered to the MNI152_T1_1mm.nii.gz image,
%       which can be found in FSL data/standard directory or in 
%       MNI_templaces within MNIparcs.
%   2. In this batch you will been three variables to describe the
%   parcellation:
%       - plabel - a short nickname for the parcellation that will be used
%       in the file naming convention within the pipeline.
%       -  pdir  - the name of the directory you were asked to create in
%       step 1 (name of the parcellation volume without the .nii.gz).
%       - pcort  - YES=1 ; NO=0; is this a cortex only parcellation? 
%       This means no cerebellum and no subcortical. Setting this to YES,
%       will result in an attempt to clean the bleeding of the parcellation 
%       into subcortical and cerebellar regions, due to transformations and 
%       dilations. 
%   3. For visualization in fMRI_B a sorting .mat file can be provided. It 
%   must contain the following variable:
%       -  ROIs  - A vector where each row represents a node and the value
%       corresponds to the grouping/network label.
%   The grouping .mat file must be in the nodal parcellation directory for
%   those nodes to be ordered according to that parcellation during the
%   visualization of the matrices.

% Initializing parcellation structure (DO NOT CHANGE)
parcs.plabel.name=struct.empty;
parcs.pdir.name=struct.empty;
parcs.pcort.true=struct.empty;
parcs.pnodal.true=struct.empty;

%% SET WHICH PARCELLATIONS YOU WANT TO USE
% shen 1_5 286 region parcellation with modified subcortical
parcs.plabel(1).name='shen_278';
parcs.pdir(1).name='shen_MNI152_org';
parcs.pcort(1).true=0;
parcs.pnodal(1).true=1;

%yeo7 resting state network parcellation
parcs.plabel(2).name='yeo7';
parcs.pdir(2).name='yeo7_MNI152';
parcs.pcort(2).true=1;
parcs.pnodal(2).true=0;

% yeo 17 resting state network parcellation
parcs.plabel(3).name='yeo17';
parcs.pdir(3).name='yeo17_MNI152';
parcs.pcort(3).true=1;
parcs.pnodal(3).true=0;
 
%%
                    %-----------------------%
                    %  SELECT GLOBAL FLAGS  %
                    %-----------------------%
                    %       =1 is ON        %
                    %       =0 is OFF       %
                    
    %  T1   %
flags.global.T1_prepare_A = 1;
flags.global.T1_prepare_B = 1;
    %  fMRI  %
flags.global.fMRI_A = 1;
flags.global.fMRI_B = 1;
    %  DWI   %
flags.global.DWI_A = 0;
flags.global.DWI_B = 0;
flags.global.DWI_C = 0;

    % Parallel %
configs.parallel = 1;
configs.UsageCPU = 0.67;

%%
                        %----------------%
                        %  T1_prepare_A  %
                        %----------------%
                        
                    %------------------------%
                    %  SELECT T1_A SUBFLAGS  %
                    %------------------------%
                    
flags.T1.dcm2niix = 1;  % dicom to nifti conversion 
    configs.T1.useCropped = 0; % use cropped field-of-view output of dcm2niix
flags.T1.denoiser = 1; % denoising
flags.T1.bet = 0; % brain extraction and mask generation (only needed for double BET)
    configs.T1.betF = 0.35;  % These are brain extraction parameters within FSL bet. (0.2)
    configs.T1.betG = 0.15; % See fsl bet help page for more details. (0.15)
flags.T1.maskCrop = 0; % Create a brain mask and use it for cropping the image (only needed for double BET)
flags.T1.anat = 1; % run FSL_anat 
    configs.T1.bias = 1; % 0 = no; 1 = weak; 2 = strong
    configs.T1.crop = 0; % 0 = no; 1 = yes (lots already done by dcm2niix)
flags.T1.bet2 = 1; % brain extraction and mask generation
    configs.T1.betF = 0.35;  % These are brain extraction parameters within FSL bet. (0.2)
    configs.T1.betG = 0.15; % See fsl bet help page for more details. (0.15)
flags.T1.re_extract = 1; % brain extraction with mask

%%
                        %----------------%
                        %  T1_prepare_B  %
                        %----------------%
      
                    %------------------------%
                    %  SELECT T1_B SUBFLAGS  %
                    %------------------------%
                    
flags.T1.reg2MNI=1;
    configs.T1.useExistingMats=0;
flags.T1.seg = 1;
    % fast segmentation smoothing (0.10 is default)
    configs.T1.segfastH = 0.25;   % was hardwired to 0.25
    % the lowest image value set to 0 when making a mask, with the
    % upper threshold (uthr), remaining fixed at 1.
    configs.T1.masklowthr = 1; % was hardwired to 1
    % flirt, dof6 cost function (default = 'corratio')
    % use 'mutualinfo' if the default fails
    configs.T1.flirtdof6cost = 'mutualinfo';  
flags.T1.parc=1;
    configs.T1.numDilReMask = 3;

%%
                          %--------------%
                          %    fMRI_A    %
                          %--------------%   
                %------------------------------------%
                %   SELECT fMRI_A SUBFLAGS (flags)   %
                %   SET fMRI_ACSFnumVoxels PARAMETERS (configs)  %
                %------------------------------------%
                
% set number of EPI sessions/scans
    configs.EPI.epiMin = 1; % minimum scan index
    configs.EPI.epiMax = 4; % maximum scan index
flags.EPI.dcm2niix = 1; % dicom import
flags.EPI.ReadHeaders = 1; % obtain pertinent scan information
    flags.EPI.UseJson = 1; % obtain pertinent scan information through json files generated by dcm2niix
flags.EPI.ICA_AROMA = 1; % ICA AROMA for motion correction, required single session melodic processed on the subjects
    flags.EPI.useExistAROMA = 1; % (optional) if ICA-AROMA has done
    flags.EPI.feat = 0; % single session melodic, required for ICA-AROMA, if done, this step can be skipped
        configs.EPI.featVersion = '3.15'; % version of feat, don't change unless the version of feat is changed
        configs.EPI.watcher = 1; % whether the featWatcher should be turn on
        configs.EPI.pre_fwhm = 6; % Melodic pre-processing: Full Width at Half Maximum of the Gaussian kernel
        configs.EPI.brainThres = 5; % Brain/background threshold, in percentage(%).
        configs.EPI.B0Unwarp = 0; % B0 field map unwarping
        configs.EPI.melodic_st = 0; % Slice timing correction
        configs.EPI.bgimage = 1; % Background image for higher-level stats overlays, don't change unless neccessary
                                    %1: Mean highres
                                    %2: First highres
                                    %3 : Mean functional
                                    %4 : First functional
                                    %5 : Standard space template
        configs.EPI.reghighres_search = 90; % Search space for registration to main structural, don't change unless neccessary
                                            %0: No search
                                            %90: Normal search
                                            %180: Full search
        configs.EPI.regstandard_search = 90; % Search space for registration to standard space, don't change unless neccessary
                                            %0: No search
                                            %90: Normal search
                                            %180: Full search 
        configs.EPI.regstandard_dof = 12; % Degrees of Freedom for registration to standard space
        configs.EPI.regstandard_nonlinear_yn = 1; % Do nonlinear registration from structural to standard space?
        configs.EPI.regstandard_nonlinear_warpres = 10; % (mm) Control nonlinear warp field resolution                            
        configs.EPI.paradigm_hp = 100; % High pass filter cutoff (s)
        configs.EPI.regstandard_res = 4; % Resampling resolution
        configs.EPI.mmthresh = 0.5; % Mixture model threshold 
    
flags.EPI.MelodicUnwarped = 0; % Use unwarped EPI obtained from melodic, only works if melodic unwarped are availabe
% if MelodicUnwarped is available, SE and Topup is unneccessary, otherwise
% it will redo unwarp agian.
flags.EPI.SpinEchoUnwarp = 0; % Requires UNWARP directory and approporiate dicoms.
flags.EPI.UseUnwarped = 0; % Use unwarped EPI if both warped and unwarped are available.
    % SPIN ECHO PAIRS (A-P, P-A) Acquistion on the Prisma
    configs.EPI.SEnumMaps = 5; % Fallback Number of PAIRS of AP and PA field maps.
    % Defaults to reading *.dcm/ima files in SE AP/PA folders
    % topup (see www.mccauslanddenter.sc.edu/cml/tools/advanced-dti - Chris Rorden's description
    % readOutTime = echoSpacing*((matrixlines4phase*partialFourier/accelrerationFactor)-1)
    % readOutTime now calculated from image data.
flags.EPI.RunTopup = 0; % 1 = Run topup (1st pass), 0 = Do not rerun if previously completed.       
    % Gradient recalled echo Field Map Acquisition
    configs.EPI.GREmagdcm = 'GREFM_MAG_DICOMS'; % MAGNITUDE Series
    configs.EPI.GREphasedcm = 'GREFM_PHASE_DICOMS'; % PHASE Series
    configs.EPI.GREbetf = 0.5; % GRE-specific bet values. Do not change
    configs.EPI.GREbetg = 0;   % GRE-specific bet input. Change if needed 
    configs.EPI.GREdespike = 1; % Perform FM despiking
    configs.EPI.GREsmooth = 3; % GRE phase map smoothing (Gaussian sigma, mm)
    configs.EPI.EPIdwell = 0.000308; % Dwell time (sec) for the EPI to be unwarped 
flags.EPI.SliceTimingCorr = 0;
    configs.EPI.UseTcustom = 1;% 1: use header-extracted times (suggested)
flags.EPI.MotionCorr = 0;
flags.EPI.RegT1 = 1;
    configs.EPI.epibetF = 0.25;
    configs.EPI.minVoxelsClust = 8; % originally hardwired to 8
flags.EPI.RegOthers = 1;
    configs.EPI.GMprobthr = 0.2;% Threshold the GM probability image
                                 % change from 0.25 to 0.2 or 0.15
flags.EPI.Mode1000 = 1;
flags.EPI.DemeanDetrend = 1;
flags.EPI.MotionRegressors = 1;
% GS is a subflag of Motion Regressors. If equals 1 global signal 
% regression is done, if 0 it is not done.
flags.EPI.GS = 0;
% Motion thresholds for scrubbing. If a hard threshold is not 
% provided, the default is set to FSL outlier guidelines of
% the 75th percentile + 1.5 times the interquartile range.
    configs.EPI.FDth = 0.20;
    configs.EPI.DVARSth = 1.28;  % Column 1 only output
% leave the below line commented unless using a value other than
% default.
%   configs.DvarCol = 1; %  DNU!!!!!!
    configs.EPI.SDth = [];
% Maximum permitted fraction of outliers to continue the analysis: 
% 0.50 - least stringent: > 50% BOLD volumes required 
% 0.40 - more conservative: > 60% BOLD good - volumes 
    configs.EPI.scrubmaxfrac = 0.50;
% Take good BOLD volumes that are adjacent (before and after) 
% to an outlier and scrub them as well! 
% 2 = tag both adjacent volumes as outliers 
% 1 = drop a good volume between two outliers (NOT yet
% implemented)
% 0 = leave "as is" (conservative)  
    configs.EPI.scrubdilate = 0;
% Exclude the first and last configs.scrubtime seconds
% from subsequent correlation calculations (FC matrix and seed)
% Suggestion: 15 (SKYRA & Standard EPI). 
% Can be longer (15-30) for Prisma/Multiband EPI (more volumes)
% DO NOT make too long to avoid losing too many BOLD volumes 
    configs.EPI.scrubtime = 15;      
flags.EPI.BandPass = 1;
    configs.EPI.fMin = .009;
    configs.EPI.fMax = .08;
flags.EPI.TissueRegressors = 1;
    configs.EPI.numCompsPCA = [0,1,3,5];% number of components to be tested (No more than 5 components, number of tests should less than 4)
flags.EPI.SpatialSmooth = 1;% The current kernel is zero, because data is reduced to ROI level.
    configs.EPI.fwhm = 0;% Full Width at Half Maximum of the Gaussian kernel
flags.EPI.ROIs = 1;
    
%% 
                    %------------------------%
                    %  MATRICES AND FIGURES  %
                    %         fMRI_B         %
                    %------------------------%
                    
                   %--------------------------%
                   %  SELECT fMRI_B SUBFLAGS  %
                   %--------------------------%
flags.EPI.FigsMotion=1; % fig 1 & 2 of motion and parameter regression
flags.EPI.FigsFC = 1;
    flags.EPI.SaveBlockMats=1; % save block-wise matrix
flags.EPI.SaveFigs = 1; % save .fig .eps .png
flags.EPI.SaveMats = 1;
                      %------------------%
                      %  SET PARAMETERS  %
                      %------------------%
configs.EPI.minVal=-0.7; % lower and upper r bound for connectome matrix visualization
configs.EPI.maxVal=0.7;
configs.EPI.numBins=40; % for the historgam number of bins
configs.EPI.Step = 1; % not sure ; something to do with scrubbing
configs.EPI.corrType = "Pearson"; % Types of correlation, 'Pearson', 'Kendall' or 'Spearman'.   
configs.EPI.minVoxelsROI = 4; % drop from 8 to 4-6 for large(r) voxels  % not used anymore 

    
%%
                        %------------------%
                        %  DWI PROCESSING  %
                        %------------------%    
    
                       %----scrubbing----------------%
                       %  SET DWI SUBFLAGS  %
                       %--------------------%
flags.DWI.dcm2niix = 1; % dicom to nifti coversion
    configs.DWI.readout = []; % if empty get from dicom; else specify value
flags.DWI.topup = 1; % FSL topup destortion field estimation
    configs.DWI.b0cut = 1; % maximum B-value to be considered B0
flags.DWI.eddy = 1; % FSL EDDY distortion correction
    configs.DWI.EDDYf = 0.3; % fsl bet threshold for b0 brain mask used by EDDY
    configs.DWI.repolON = 1; % use eddy_repol to interpolate missing/outlier data
flags.DWI.DTIfit = 1; % Tensor estimation and generation of scalar maps
    configs.DWI.DTIfitf = 0.4; % brain extraction (FSL bet -f) parameter    
    
%%
                    %---------------------------%
                    %  STRUCTURAL CONNECTIVITY  %
                    %---------------------------%  
                    
                      %----------------------%
                      %  SET sCONN SUBFLAGS  %
                      %----------------------%
flags.DWI.reg2T1 = 1; % register diffusion data to T1 (dof6, bbr)
flags.DWI.tissueMasks = 1; % create masks 4 seeds, fibers, etc
flags.DWI.Camino = 1; % run camino processing
    configs.DWI.CaminoReset = 1; % erase existing camino directory and start over
    configs.DWI.HeapSizeCamino = 16384; % RAM allocation for camino
flags.DWI.getTensor = 1; % get tensor and fiber orientation data
    configs.DWI.order2Threshold = 12;
    configs.DWI.order4Threshold = 6;% Decrease value get more two tensor voxels
    configs.DWI.order6Threshold = 0;% Do not change this unless you are going beyond two components in deterministic modeling.
    configs.DWI.clusterMinSizeMultiTensor = 8;% Minimum cluster size so that voxels are modeled with multi-tensor. They go back to single-tensor otherwise
flags.DWI.Deterministic = 1; % deterministic whole brain tractography   
    configs.DWI.seedsWMborder = 1; % use GM WM interface for seeds
    configs.DWI.CURVthresh = 45; % maximum accepted turn ange
    configs.DWI.stepSize = 1; % tracking step, relative to voxels
flags.DWI.FiberFiles = 1; % generate trk and vtk files for visualization
    configs.DWI.LengthMin = 8; % minimum accepted fiber length
    configs.DWI.LengthMax = 180; % maximum accepted fiber length
    configs.DWI.FAthresh = 0.1; % minimum FA to allow tracking
flags.DWI.genMats = 0; % generate connectivity matrices

%%
                    %------------------------------%
                    %  PROBABILISTIC TRACTOGRAPHY  %
                    %------------------------------%  
                    
                       %--------------------------%
                       %  SET prob_conn SUBFLAGS  %
                       %--------------------------%
flags.DWI.bedpostXprep = 0; % create directory and links for bedpostX input
    configs.DWI.outsourceBpX = 0; % copy inputs to group directory for processing elsewhere (Karst)
    paths.DWI.groupDir =''; % name of group directory (use absolute path)
flags.DWI.runBedpostX = 0; % this takes ~10-20hrs per subject
    params.DWI.fibres = 2; % max number of crossing fibres per voxel
    params.DWI.weight = 1;
    params.DWI.burnin = 1000;
flags.DWI.reg2DWI = 0;
    configs.DWI.Warp = 0; % if on dof12 is done after 6 (good if no distortion field available)
flags.DWI.probtrackX = 1; % requires T1 and reg2DWI to have been done.
    configs.DWI.importBpX = 0;
    params.DWI.SeedType = [2 2]; % [seed target] types. If target is blank (e.g. [1]) tracking is unconstrained.
                                 % seed type is required.
              % 1 - use a node from parcellation as input 
                    configs.DWI.SeedIdx = []; % vector comma separated
                    configs.DWI.TargetIdx = [];
              % 2 - use an MNI space seed mask (1mm 152MNI space)
                  % The following can be 3D binary, 3D indexed, or 4D where 4th dimention is ROIs.
                    paths.DWI.MNIseedImage = '/datay1/datadir/seed/harvard-frontal_combined.nii.gz';
                    paths.DWI.MNItargetImage = '/datay1/datadir/L_striatum.nii.gz';
              % Labels (optional) - Order of labels (column format) in text file must be the
                  % same as image/seed order/index, e.g. 1=first label name; 2=second label name, etc.
                  % Labels must be a single continous string per line (no spaces)
                    paths.DWI.seedLabels = '/datay1/datadir/seed/harvard-frontal-labels.txt';
                    paths.DWI.targetLabels = '';
    params.DWI.network = 1; %BE CAREFUL!
        configs.DWI.netName = 'frontostriatal';
        % Turning on the network option with compine seed and target
        % options into a single set of seeds. The output of probtrackX2
        % with then be a tract that is the combination of all seeds. An
        % additional text file will be output with a streamline count
        % matrix among all pairs of seeds.
    configs.DWI.curvature = 0.2;
    configs.DWI.maxSteps = 2000;
    configs.DWI.stepLength = 0.5;
    configs.DWI.numSamples = 5000;
    configs.DWI.disttresh = 0; %mm
    configs.DWI.sampvox = 0.0;