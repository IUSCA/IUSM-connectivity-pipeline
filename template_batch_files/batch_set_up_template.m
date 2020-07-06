%                               BATCH_SET_UP
% This contains all flags and configs required by the pipeline to process
% the data. You may edit this as necessary depending on what portions of
%                       the pipeline you wish to run.
%
%      Evgeny Chumin, Indiana University, Bloomington, 2019
%      Zikai Lin, Indiana University School of Medicine, 2018
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
% parcs.plabel.name=struct.empty;
% parcs.pdir.name=struct.empty;
% parcs.pcort.true=struct.empty;
% parcs.pnodal.true=struct.empty;

%% SET WHICH PARCELLATIONS YOU WANT TO USE
% Schaefer parcellation of yeo17 into 200 nodes
% parcs.plabel(1).name='schaefer200_yeo17';
% parcs.pdir(1).name='Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm';
% parcs.pcort(1).true=1;
% parcs.pnodal(1).true=1;
% 
% % Schaefer parcellation of yeo17 into 300 nodes
% parcs.plabel(2).name='schaefer300_yeo17';
% parcs.pdir(2).name='Schaefer2018_300Parcels_17Networks_order_FSLMNI152_1mm';
% parcs.pcort(2).true=1;
% parcs.pnodal(2).true=1;
% 
% % DKI parcellation
% parcs.plabel(3).name='DKIcort';
% parcs.pdir(3).name='DKIcort';
% parcs.pcort(3).true=1;
% parcs.pnodal(3).true=1;

parcs.plabel(1).name='shen';
parcs.pdir(1).name='blank';
parcs.pcort(1).true=0;
parcs.pnodal(1).true=1;
 
%%
                    %-----------------------%
                    %  SELECT GLOBAL FLAGS  %
                    %-----------------------%
                    %       =1 is ON        %
                    %       =0 is OFF       %
                    
    %  T1   %
flags.global.T1_prepare_A = 0;
flags.global.T1_prepare_B = 0;
    %  fMRI  %
flags.global.fMRI_A = 0;
    %  DWI   %
flags.global.DWI_A = 0;
flags.global.DWI_B = 1; % not fully tested; run at your own risk
flags.global.DWI_C = 0; % depricated; may not work

%%
                        %----------------%
                        %  T1_prepare_A  %
                        %----------------%
                    %------------------------%
                    %  SELECT T1_A SUBFLAGS  %
                    %------------------------%
                    
flags.T1.dcm2niix = 0;  % dicom to nifti conversion 
    configs.T1.useCropped = 0; % use cropped field-of-view output of dcm2niix
flags.T1.denoiser = 0; % denoising % 2 - do not denoise, copy input to output
                                   % 1 - denoise MRI
                                   % 0 - skip if already ran(do not copy input to output)
flags.T1.anat = 1; % run FSL_anat 
    configs.T1.bias = 1; % 0 = no; 1 = weak; 2 = strong;
    configs.T1.crop = 1; % 0 = no; 1 = yes (lots already done by dcm2niix if used)
flags.T1.bet = 0; % brain extraction and mask generation
    configs.T1.antsTemplate = 'NKI'; % 'MICCAI' or 'NKI' or 'IXI' or 'bet'
    configs.T1.betF = .2;  % These are brain extraction parameters within FSL bet. (0.2)
    configs.T1.betG = -0.02; % See fsl bet help page for more details. (0.15)
% currently testing ANTS, which does not require bet inputs
flags.T1.re_extract = 1; % brain extraction with mask

%%
                        %----------------%
                        %  T1_prepare_B  %
                        %----------------%
                    %------------------------%
                    %  SELECT T1_B SUBFLAGS  %
                    %------------------------%
                    
flags.T1.reg2MNI=1;
    configs.T1.useExistingMats=1;
    configs.T1.useMNIbrain=1; % use MNI152T1_brain (=1) rather than MNI152T1 (=0)
    configs.T1.fnirtSubSamp='4,2,2,1'; %subsampling for fnirt registration to MNI; fnirt default string = '4,2,1,1' pipeline default = '4,2,2,1'
flags.T1.seg = 1;
    % fast segmentation smoothing (0.10 is default)
    configs.T1.segfastH = 0.25;   % was hardwired to 0.25
    % flirt, dof6 cost function (default = 'corratio')
    % use 'mutualinfo' if the default fails
    configs.T1.flirtdof6cost = 'mutualinfo'; % 'corratio'-fsl default; 'mutualinfo'-recommended 
flags.T1.parc=1;
    configs.T1.numDilReMask = 2;
    configs.T1.addsubcort=1; % add FSL subcortical to cortical parcellations ONLY
                             % to nodal parcellations individuals regions
                             % are added, while to others they are added as
                             % a single subcortical network.
%%
                          %--------------%
                          %    fMRI_A    %
                          %--------------%   
                     %---------------------------%
                     %   SELECT fMRI_A SUBFLAGS  %
                     %---------------------------%
                            % DATA IMPORT
% Number of EPI sessions is determined by number of EPI directories. If you
% wish to only process select sessions, set them via epiMin and epiMax.
    configs.EPI.epiMin = 1; % minimum scan index
    configs.EPI.epiMax = 1; % maximum scan index
flags.EPI.dcm2niix = 0; % dicom import
flags.EPI.ReadHeaders = 1; % obtain pertinent scan information
    configs.EPI.UseJson = 1; % obtain pertinent scan information through json files generated by dcm2niix
%-------------------------------------------------------------------------%
                            % DISTORTION CORRECTION
                %--------------------------------------------%
                %  Spin Echo Field Map Acquisition
flags.EPI.SpinEchoUnwarp = 0; % Requires UNWARP directory and approporiate dicoms.
    % SPIN ECHO PAIRS (A-P, P-A) Acquistion on the Prisma
    configs.EPI.SEnumMaps = 3; % Fallback number of PAIRS of AP and PA field maps.
    % SEnumMaps is delaulted to if fslinfo cannot determine number of
    % images.
                %--------------------------------------------%
                % Gradient recalled echo Field Map Acquisition
    configs.EPI.GREbetf = 0.5; % GRE-specific bet values. Do not change
    configs.EPI.GREbetg = 0;   % GRE-specific bet input. Change if needed 
    configs.EPI.GREdespike = 1; % Perform FM despiking
    configs.EPI.GREsmooth = 3; % GRE phase map smoothing (Gaussian sigma, mm)
    configs.EPI.EPIdwell = 0.000308; % Dwell time (sec) for the EPI to be unwarped 
%-------------------------------------------------------------------------%
                    % PREPROCESSING & ANATOMY REGISTRATION
flags.EPI.SliceTimingCorr = 1; % recommended for TR > 1.2s or non multiband data
    configs.EPI.minTR = 1.7; % skip if TR below minTR; good for multiple acquisition datasets
    configs.EPI.UseTcustom = 1;% 1: use header-extracted times (suggested)
flags.EPI.MotionCorr = 1;
flags.EPI.RegT1 = 1;
    configs.EPI.epibetF = 0.3;
    configs.EPI.minVoxelsClust = 8; % originally hardwired to 8
flags.EPI.RegOthers = 1;
    configs.EPI.GMprobthr = 0.2;% Threshold the GM probability image
                                 % change from 0.25 to 0.2 or 0.15
flags.EPI.IntNorm4D = 1; % Intensity normalization to global 4D mean of 1000
%-------------------------------------------------------------------------%
                       % MOTION AND OUTLIER CORRECTION
%-------------------------------------------------------------------------% 
% This section must be ran as a whole!
flags.EPI.NuisanceReg = 2;
    % 0 - skip
    % 1 - ICA-based denoising; WARNING: This will smooth your data.
    % 2 - Head Motion Parameter Regression
        configs.EPI.numReg = 24; % 12 (orig and deriv) or 24 (+ sq of 12)
    flags.EPI.PhysReg = 2; %physiological regressors
        % 0 - skip
        % 1 - aCompCorr; PCA based CSF and WM signal regression (up to 5
        %     components)
            configs.EPI.numPC = 5; % 1-5; the maximum and recommended number is 5 % leave empty to putput all
        % 2 - mean WM and CSF signal regression
            configs.EPI.numPhys = 8; % 2-orig; 4-orig+deriv; 8-orig+deriv+sq
    flags.EPI.GS = 1; % global signal regression 
        configs.EPI.numGS = 4; % 1-orig; 2-orig+deriv; 4-orig+deriv+sq        
    flags.EPI.DemeanDetrend = 1;
    flags.EPI.BandPass = 1;
        configs.EPI.fMin = .008; %Parkes 2018 & Satterthwaite 2013 NeuroImage
        configs.EPI.fMax = .08;   
    flags.EPI.scrubbing = 0; % perform scrubbing based on FD and DVARS criteria
        configs.EPI.FDcut = []; % [] uses fsl outlier critera
        configs.EPI.DVARScut = [];
    flags.EPI.ROIs = 1; % extract parcellation timeseries
%-------------------------------------------------------------------------%   
%%
                        %------------------%
                        %  DWI PROCESSING  %
                        %------------------%    
    
                       %--------------------%
                       %  SET DWI SUBFLAGS  %
                       %--------------------%
flags.DWI.dcm2niix = 0; % dicom to nifti coversion
    configs.DWI.readout = []; % if empty get from dicom; else specify value
flags.DWI.topup = 0; % FSL topup destortion field estimation
    configs.DWI.b0cut = 1; % maximum B-value to be considered B0
flags.DWI.eddyPREP = 0; % FSL EDDY distortion correction
    configs.DWI.EDDYf = 0.3; % fsl bet threshold for b0 brain mask used by EDDY
flags.DWI.eddyRUN = 0;   
    configs.DWI.repolON = 1; % use eddy_repol to interpolate missing/outlier data
flags.DWI.DTIfit = 0; % Tensor estimation and generation of scalar maps
    configs.DWI.DTIfitf = 0.3; % brain extraction (FSL bet -f) parameter    
    
%%
                    %---------------------------%
                    %  STRUCTURAL CONNECTIVITY  %
                    %---------------------------%  
                    
                      %----------------------%
                      %  SET sCONN SUBFLAGS  %
                      %----------------------%
flags.DWI.regT1_2DWI = 0;
flags.DWI.regParc = 1;
flags.DWI.MRtrix = 0;
    configs.DWI.seeding = 'wm'; % wm - seed white matter mask
                                % dyn - dynamic random seeding
    configs.DWI.numSeeds_perWMvoxel = 2;
flags.DWI.connMatrix = 1; % generate connectivity matrices

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
flags.DWI.probtrackX = 0; % requires T1 and reg2DWI to have been done.
    configs.DWI.importBpX = 0;
    params.DWI.SeedType = [2 2]; % [seed target] types. If target is blank (e.g. [1]) tracking is unconstrained.
                                 % seed type is required.
              % 1 - use a node from parcellation as input 
                    configs.DWI.SeedIdx = []; % vector comma separated
                    configs.DWI.TargetIdx = [];
              % 2 - use an MNI space seed mask (1mm 152MNI space)
                  % The following can be 3D binary, 3D indexed, or 4D where 4th dimention is ROIs.
                    paths.DWI.MNIseedImage = '/path/to/your/MNI_seeds.nii.gz';
                    paths.DWI.MNItargetImage = '/path/to/your/MNI_targets.nii.gz';
              % Labels (optional) - Order of labels (column format) in text file must be the
                  % same as image/seed order/index, e.g. 1=first label name; 2=second label name, etc.
                  % Labels must be a single continous string per line (no spaces)
                    paths.DWI.seedLabels = '/path/to/your/seed-labels.txt';
                    paths.DWI.targetLabels = '/path/to/your/target-labels.txt';
    params.DWI.network = 1; %BE CAREFUL!
        configs.DWI.netName = 'name-your-network';
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