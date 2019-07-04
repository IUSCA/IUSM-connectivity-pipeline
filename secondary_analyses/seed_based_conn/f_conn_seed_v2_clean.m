clearvars
clc
%% ------------------------------- PATHS ------------------------------- %%

% Exporting FSL path
FSLsetup = 'FSLDIR=/usr/local/fsl; . ${FSLDIR}/etc/fslconf/fsl.sh; PATH=${FSLDIR}/bin:${PATH}; export FSLDIR PATH';
[~,result]=system(FSLsetup);
% Addings connectome spripts path (Includes NIfTI toolbox)

addpath(genpath('/data03/CONNECTIVITY/connectivity_pipeline_v1.7/connectome_scripts_v1.7'))

%-------------------------------------------------------------------------%

% Creating a directory for results
paths.results = fullfile(pwd,'..','seed_results');
    if ~exist(paths.results,'dir')
        mkdir(paths.results)
    else
        warning('Results will be written to existing directory:\n %s',paths.results)
    end
    
%% --------------------- USER DEFINED INPUT DATA ----------------------- %%

% Study ID / prepend to the subject IDs
StudyID = 'BRACS';

% Data directory (contains subject directories)
paths.data = '/data03/CONNECTIVITY/connectivity_pipeline_v1.7/datadir';

% Name of directories containining EPI data for each subject
ScanID = 'EPI2'; % 'REST' or 'EPI1' or 'whatever'
GSregID = 'GSreg_yes'; % 'GSreg_yes' or 'GSreg_no' (GS regression on/off)
PCAID = 'PCA3';  % 'PCA0', 'PCA1', 'PCA3', or 'PCA5' (see connectome scripts)
paths.resultscp = fullfile(paths.results,strcat(ScanID,'_',GSregID,'_',PCAID));
if ~exist(paths.resultscp,'dir')
    mkdir(paths.resultscp)
end

% List of subject numbers

slist = [27; 186; 228; 341; 363; 914; 1063; 1069; 1079; 1106; 1139; ...
        1156; 1167; 1172; 1186; 1313; 1331; 1342; 1349; 1364; 1369;...
        1372; 1377; 1387; 1397; 1401; 1412; 1422; 1473; 1480; 1482;...
        1488; 1495; 1496; 1499; 1504; 1505; 1512; 1520; 1522; 1535;...
        1536; 1551; 1575; 1865; 1871; 1887; 1888; 1896; 1897; 1899;...
        1906; 1907; 1913; 1922; 1927; 1928; 1946; 1947; 1951; 1966];
% slist = [ 27];
% if EPI-to-MNI transformation matrices already exist (previous seed), 
% should they be recalculated for subsequent seeds?
% 1 = Yes (slows down the code but helpful if starting from scratch)
% 0 = No (calculate once, then use those matrices for subsequent seeds)
configs.TransMatOverwrite = 0;

%% --------------------------- DEFINE SEEDS ---------------------------- %%

% Provide Shen IDs; Indicate parcellation version in CONFIGS
% list.regions = [25; 74; 81; 92; 114; 134; 151];
%  list.regions = [114];
list.regions = [];

% Provide MNI coordinates X Y Z (tab delimited) 3-column test file as origins 
% for generation of seeds.
seedCoordList = ...
        '/data03/CONNECTIVITY/connectivity_pipeline_v1.7/secondary_analyses/seed_based_setup/testing_seeds.txt';
    %'/localscratch/jenya_connectome_scripts/datadir/testing_seeds.txt';

% Define Seed parameters:
% Shape     % =0 Sphere
            % =1 Rectangle
configs.seedShape =1; 

% Define Seed Size
if configs.seedShape == 0
    configs.seedRad = 5;% Radius in mm of the sphere
else
    %           Must be odd values!
    configs.seedBox = [5 5 5]; %[ X   Y   Z ] dimentions in voxels of the box. 
    %                            [L-R A-P S-I]
end

%Provide a path to an image you wish to use as a seed, as well as a label
%for the seed.
    % This image will binarized in the code, therefore if you wish to
    % impose a threshold, do so before providing image as input.
    % The input must be in MNI152_T1 space.
    seedLabel = 'L_NAcc';
    seedImage = ...
    '/data03/CONNECTIVITY/connectivity_pipeline_v1.7/secondary_analyses/seed_based_setup/L_NAcc_2mm.nii.gz';

%% ----------------------------- CONFIGS ------------------------------- %%

% Use only Gray Matter voxels 
configs.GMmask = 0; % 1 = "On"; 0 = "Off" 

% Number of volumes dropped off the start and end of EPI series.
configs.dropTR = 15;

% this is the level of erosion on Shen-seeds. 1 does nothing. 
% 3 most likely recommended. 5 might produce empty seeds. 
% even numbers will produce uneven erosions, but are possible.
configs.erodeLevel = 3;

% Specify the resolution of MNI space seed conn results: 1mm or 2mm.
configs.MNIres = 2;

% Smoothing kernel Full-width Half-max (mm)
configs.fwhm = 6;

%% ------------------- CHOOSE WHICH INPUT TO USE ----------------------- %%
    % Shen region IDs as input
                    flags.runShens = 0;
                    
    % List of MNI coordinates as input
                    flags.runMNIcoords = 0;
                    
    % MNI space image as input                
                    flags.runImage = 1; % Currently being tested
%%
                    
                    
                    
                    
                    
                    
                    
                    
%% --------------ALL CODE BELOW THIS LINE IS NO TOUCHIE!---------------- %%

% Generate a list of subjects based on provided list and existing data
% directories.
for i=1:size(slist,1)
    StudyNo = sprintf('%0.4d',slist(i));
    studyname = strcat(StudyID,StudyNo);
    if exist(fullfile(paths.data,studyname),'dir')
        list.subjects(i) = struct('name',studyname);
    end
end
%-------------------------------------------------------------------------%

% ---------------------Looping Across Subjects--------------------------- %
for i=1:length(list.subjects)
    list.subjects(i).name
    % Set subject specific paths
    paths.EPI = fullfile(paths.data,list.subjects(i).name,ScanID);
    paths.GSPCA = fullfile(paths.EPI,GSregID,PCAID);
    paths.seed = fullfile(paths.GSPCA,'SeedConnectivity');
    if ~exist(paths.seed,'dir')
        mkdir(paths.seed)
    else
        warning('Results for %s will be written to existing directory:\n %s',...
            list.subjects(i).name,paths.seed)
    end
    
    % Read in the EPI 4D volume for the subject
    subjNii = MRIread(fullfile(paths.GSPCA,'9_epi.nii.gz'));
    % Drop early and late volumes
    subjVOL = subjNii.vol(:,:,:,configs.dropTR+1:end-configs.dropTR);clear subjNii;
    % Calculate numver of remaining time points
    numVols = size(subjVOL,4);
    
    % Read in the EPI space parcellation
    subjPARC = MRIread(fullfile(paths.EPI,'rT1_GM_parc_shen_clean.nii.gz'));
    % Create mask from parcellation
    subjMaskGM = subjPARC.vol > 0;
    % Create a Temp volume structure for writing out NIfTI files
    Temp = subjPARC;
    
% -----------Prepare transformation matrices for MNI space--------------- %
    
% EPI -> T1 linear dof6 and bbr (inverse of T1 -> EPI)
    fileMat = fullfile(paths.EPI,'T1_2_epi_dof6_bbr.mat');
    fileOut = 'epi_dof6_bbr_2_T1.mat';
    fileMatInv = fullfile(paths.EPI,fileOut);
    if ~exist(fileMatInv,'file') || configs.TransMatOverwrite == 1;
        sentence = sprintf('convert_xfm -omat %s -inverse %s',fileMatInv,fileMat);
        [status,result] = system(sentence);
    else
        sprintf('Using already available %s:\n',fileOut)
    end

% Combine with T12MNI_dof6
    paths.reg = fullfile(paths.data,list.subjects(i).name,'T1/registration');
    fileMat1 = fullfile(paths.reg,'T12MNI_dof6.mat');
    fileMat2 = fullfile(paths.EPI,'epi_dof6_bbr_2_T1.mat');
    fileOut = 'epi_2_MNI_dof6.mat';
    fileMatJoint = fullfile(paths.EPI,fileOut); 
    if ~exist(fileMatJoint,'file') || configs.TransMatOverwrite == 1;
        sentence = sprintf('convert_xfm -omat %s -concat %s %s',fileMatJoint,fileMat1,fileMat2);
        [status,result] = system(sentence);
    else
        sprintf('Using already available %s:\n',fileOut)
    end
    
% Finally, apply T12MNI_dof12 and make an inverse
    fileMat1 = fullfile(paths.reg,'T12MNI_dof12.mat');
    fileMat2 = fullfile(paths.EPI,'epi_2_MNI_dof6.mat');
    fileOut = 'epi_2_MNI_final.mat';
    fileMatJoint = fullfile(paths.EPI,fileOut);
    if ~exist(fileMatJoint,'file') || configs.TransMatOverwrite == 1;
        sentence = sprintf('convert_xfm -omat %s -concat %s %s',fileMatJoint,fileMat1,fileMat2);
        [status,result] = system(sentence);
    else
        sprintf('Using already available %s:\n',fileOut)
    end

    fileMat = fullfile(paths.EPI,'epi_2_MNI_final.mat');
    fileOut = 'MNI_2_epi_final.mat';
    fileMatInv = fullfile(paths.EPI,fileOut);
    if ~exist(fileMatInv,'file') || configs.TransMatOverwrite == 1;
        sentence = sprintf('convert_xfm -omat %s -inverse %s',fileMatInv,fileMat);
        [status,result] = system(sentence);
    else
        sprintf('Using already available %s:\n',fileOut)
    end
    % T1 -> MNI nonlinear warp  - generated in the pipeline.
%-------------------------------------------------------------------------%    

%% ------------------Looping Across Shen Seeds-------------------------- %%
if flags.runShens == 1
    for j=1:length(list.regions)
        list.regions(j)
        % Isolate seed from parcellation
        seedMask = subjPARC.vol == list.regions(j);
        % Erode the seed.
        seedMask = imerode(seedMask,ones(configs.erodeLevel,configs.erodeLevel,configs.erodeLevel));
        % If the seed is gone, repeat with one less erosion.
        if nnz(seedMask)<2
            seedMask = subjPARC.vol == list.regions(j);
            seedMask = imerode(seedMask,ones(configs.erodeLevel-1,configs.erodeLevel-1,configs.erodeLevel-1));
        end
        % If the seed is still gone, skip erosion.
        if nnz(seedMask)<2
            seedMask = subjPARC.vol == list.regions(j);
        end
        
        % Reshape 4D mat to a 2D, where rows=voxels and colums=volumes/time
        voxels = reshape(subjVOL,[size(subjVOL,1)*size(subjVOL,2)*...
            size(subjVOL,3),numVols]);
        % Extract time-series of seed voxels.
        seed_voxels=voxels(seedMask,:);
        % Get mean time-series of the seed
        seed=mean(seed_voxels,1);
        
        % Correlate mean seed time-series with whole-brain
        C.vector = corr(seed',voxels');
        % Isolate mean seed time-series correlations with seed voxels.
        S.vector = C.vector;
        S.vector(~seedMask) = 0.;
        
        % Fisher's Z transformation of correlations
        C.vector = 0.5*log((1+C.vector)./(1-C.vector));
        S.vector = 0.5*log((1+S.vector)./(1-S.vector));
        
% -----------Write NIfTI volumes Fisher's Z of correlations-------------- %

             % ----------------Whole-Brain ---------------%
        C.vol = nan(size(subjMaskGM));
        C.vol = reshape(C.vector,size(C.vol));
        Temp.vol = C.vol;
        MRIwrite(Temp,fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_wholebrain_epi.nii.gz',list.regions(j))));
        
             % ------------------Seed-Only ---------------%
        S.vol = zeros(size(subjMaskGM));
        S.vol = reshape(S.vector,size(S.vol));
        Temp.vol = S.vol;
        MRIwrite(Temp,fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_seedonly_epi.nii.gz',list.regions(j))));
        
% --------------------Transform images into MNI space-------------------- %          

             % ----------------Whole-Brain ---------------%
        % Input Fisher's Z image in EPI space.
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_wholebrain_epi.nii.gz',list.regions(j)));
        % Output MNI space image name.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_wholebrain_MNI.nii.gz',list.regions(j)));
        % Transformation field and matrix to MNI space.
        fileMat = fullfile(paths.EPI,'epi_2_MNI_final.mat');
        fileWarp = fullfile(paths.reg,'T12MNI_warp.nii.gz');
        % Select reference image based on specified resolution.
        if configs.MNIres == 1
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
        elseif configs.MNIres == 2
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
        else
            error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
        end
        % Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        
             % ------------------Seed-Only ---------------%
        % Input seed Fisher's Z
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_seedonly_epi.nii.gz',list.regions(j)));
        % Output seed MNI space image.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_seedonly_MNI.nii.gz',list.regions(j)));
        % Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);    
        
% --------------Copy images into group results directory----------------- %

             % ----------------Whole-Brain ---------------%
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_wholebrain_MNI.nii.gz',list.regions(j))); 
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed%03d_eroded_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,list.regions(j)));
        sentence = sprintf('cp %s %s',fileIn,fileOut);
        [~,result]=system(sentence);
        
             % ------------------Seed-Only ---------------%
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed%03d_eroded_seedonly_MNI.nii.gz',list.regions(j)));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed%03d_eroded_seedonly_MNI.nii.gz',...
        list.subjects(i).name,list.regions(j)));
        sentence = sprintf('cp %s %s',fileIn,fileOut);
        [~,result]=system(sentence);
        
% ----------------------Smooth the MNI space data------------------------ %

                % --------convert FWHM to sigma-------- %
        configs.sigma = (configs.fwhm/sqrt(8*log(2)));
        
             % ----------------Whole-Brain ---------------%
        % Read in MNI space Fisher's Z image
        brainMNI = MRIread(fullfile(paths.resultscp,...
            sprintf('%s_conn_seed%03d_eroded_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,list.regions(j))));
%         % Replace NaN with zeros
        brainMNI.vol(isnan(brainMNI.vol))=0;
%         % Overwrite the NaN containing image with zero containing image
        MRIwrite(brainMNI,fullfile(paths.resultscp,...
            sprintf('%s_conn_seed%03d_eroded_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,list.regions(j))));
%         % Smooth the data
        fileIn = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed%03d_eroded_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,list.regions(j)));
        fileSmoothed = fullfile(paths.resultscp,...
            sprintf('%.0ffwhm_%s_conn_seed%03d_eroded_wholebrain_MNI.nii.gz',...
        configs.fwhm,list.subjects(i).name,list.regions(j)));
        sentence = sprintf('fslmaths %s -s %.8f %s',fileIn,configs.sigma,fileSmoothed);
            [~,result]=system(sentence);
    end
%-------------------------------------------------------------------------%
end
%% ----------Transforming Seed MNI Coordiantes into EPI Space----------- %%
if flags.runMNIcoords == 1
    
    % Read in MNI seed coordinate list from input into a matrix.
    if ~exist(seedCoordList,'file')
        error('No text file with MNI coordinates provided. Terminating...')
    else
        MNIcoords = dlmread(seedCoordList,'\t');
    end
    
    % ----------------------MNI -> T1 warp--------------------- %
    fileSource = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz'; % used by the pipeline - DO NOT CHANGE
    fileDest = fullfile(paths.data,list.subjects(i).name,...
        'T1/registration/T1_dof12.nii.gz');
    filewarp = fullfile(paths.data,list.subjects(i).name,...
        'T1/registration/MNI2T1_warp.nii.gz');
    fileTmp1 = fullfile(paths.seed,'seed_coords_unwarped.tmp');
    fileOut1 = fullfile(paths.seed,'seed_coords_unwarped.txt');
    sentence = sprintf('cat %s | img2imgcoord -src %s -dest %s -warp %s -mm > %s',...
        seedCoordList,fileSource,fileDest,filewarp,fileTmp1);
    [~,result]=system(sentence);
% grep -v "Coordinates in Destination volume (in mm)" /data04/CONNECTIVITY/connectivity_pipeline_v1.5/datadir/BRACS0341/EPI2/GSreg_yes/PCA5/SeedConnectivity/seed_coords_unwarped.txt|awk '{print $1","$2","$3}'
    str1 = 'grep -v "Coordinates in Destination volume (in mm)"';
    str2 = '|awk ''{print $1"\t"$2"\t"$3}''';
    test1 = sprintf('%s %s %s', str1, fileTmp1, str2);
    sentence = sprintf('%s > %s',test1,fileOut1); 
    [~,result]=system(sentence);
       
    % --------------T1 -> EPI dof6 and bbr------------------ %
    fileSource = fullfile(paths.data,list.subjects(i).name,...
        'T1/registration/T1_dof12.nii.gz');
%     fileDest = fullfile(paths.EPI,'rT1_brain_dof6.nii.gz');
    fileDest = fullfile(paths.EPI,'2_epi_meanvol.nii.gz');
    filedof6bbr = fullfile(paths.EPI,'MNI_2_epi_final.mat');
    fileTmp4 = fullfile(paths.seed,'seed_coords_EPIspace.tmp');    
    fileOut4 = fullfile(paths.seed,'seed_coords_EPIspace.txt');
    sentence = sprintf('cat %s | img2imgcoord -src %s -dest %s -xfm %s -mm > %s',...
        fileOut1,fileSource,fileDest,filedof6bbr,fileTmp4);
    [~,result]=system(sentence); 
    str1 = 'grep -v "Coordinates in Destination volume (in mm)"';
    str2 = '|awk ''{print $1"\t"$2"\t"$3}''';
    test1 = sprintf('%s %s %s', str1, fileTmp4, str2);
    sentence = sprintf('%s > %s',test1,fileOut4); 
    [~,result]=system(sentence);   
    
    % Convert obtained EPI space seed coordinates from mm to voxel
    fileIMG = fullfile(paths.EPI,'rT1_brain_mask.nii.gz');
    fileOut5 = fullfile(paths.seed,'seed_voxelCoords_EPIspace.txt');
    sentence = sprintf('cat %s | std2imgcoord -img %s -std %s -vox > %s',...
        fileOut4,fileIMG,fileIMG,fileOut5);
    [~,result]=system(sentence);
    
%% -----------Generate point ROIs for each seed coordinate-------------- %%
% ------then grow them to the shape and size specified in configs-------- %

    % ---------------------Generating points-------------------- %
    
    % Read in voxel EPI space coordinates
    ImageIN = fullfile(paths.EPI,'rT1_brain_mask.nii.gz');
    voxelCoords = dlmread(fileOut5);
    % For each set of coordinates:
    for l=1:size(voxelCoords,1)
        % Generate a single point at the coordinate
        SeedName = sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_pointROI.nii.gz',...
            MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3));
        SeedOut = fullfile(paths.seed,SeedName);
        sentence = sprintf('fslmaths %s -roi %s 1 %s 1 %s 1 0 1 %s'...
            ,ImageIN...
            ,num2str(voxelCoords(l,1)),num2str(voxelCoords(l,2)),num2str(voxelCoords(l,3))...
            ,SeedOut);
        [~,result]=system(sentence);
        
    % ---------------------Growing the seeds-------------------- %
        if configs.seedShape == 0 % Sphere
            SeedName2 = sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphereROI.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3));
            SeedOut2 = fullfile(paths.seed,SeedName2);
            sentence = sprintf('fslmaths %s -kernel sphere %d -fmean -bin %s',...
                SeedOut,configs.seedRad,SeedOut2);
            [~,result]=system(sentence);
            
        else % Rectangle
            SeedName2 = sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboidROI.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3));
            SeedOut2 = fullfile(paths.seed,SeedName2);
            sentence = sprintf('fslmaths %s -kernel boxv3 %s %s %s -fmean -bin %s',...
                SeedOut,...
                num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)),...
                SeedOut2);
            [~,result]=system(sentence);
        end
        
    % Set seedMask as the newly generated ROI
    seedMask = MRIread(SeedOut2);
    if configs.GMmask == 1  % Use only GM ROI voxels
        seedMask = logical(seedMask.vol.*subjMaskGM); 
    else                    % Include ALL ROI voxels
        seedMask = logical(seedMask.vol);
    end
    
    % Reshape 4D mat to a 2D, where rows=voxles and colums=volumes/time
    voxels = reshape(subjVOL,[size(subjVOL,1)*size(subjVOL,2)*...
        size(subjVOL,3),numVols]);
    % Extract time-series of seed voxels.
    seed_voxels=voxels(seedMask,:);
    % Get mean time-series of the seed
    seed=mean(seed_voxels,1);
    
    % Correlate mean seed time-series with whole-brain
    C.vector = corr(seed',voxels');
    % Isolate mean seed time-series correlations with seed voxels.
    S.vector = C.vector;
    S.vector(~seedMask) = 0.;
    
    % Fisher's Z transformation of correlations
    C.vector = 0.5*log((1+C.vector)./(1-C.vector));
    S.vector = 0.5*log((1+S.vector)./(1-S.vector));

% -----------Write NIfTI volumes Fisher's Z of correlations-------------- %

             % ----------------Whole-Brain ---------------%
    C.vol = nan(size(subjMaskGM));
    C.vol = reshape(C.vector,size(C.vol));
    Temp.vol = C.vol;
    
    if configs.seedShape == 0 % Sphere
    MRIwrite(Temp,fullfile(paths.seed,...
        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_epi.nii.gz',...
            MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
            configs.seedRad)));
    else % Cuboid
    MRIwrite(Temp,fullfile(paths.seed,...
        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_epi.nii.gz',...
            MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
            num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
    end
    
             % ------------------Seed-Only ---------------%
    S.vol = zeros(size(subjMaskGM));
    S.vol = reshape(S.vector,size(S.vol));
    Temp.vol = S.vol;
    
    if configs.seedShape == 0 % Sphere
    MRIwrite(Temp,fullfile(paths.seed,...
        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_epi.nii.gz',...
            MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
            configs.seedRad)));
    else % Cuboid
    MRIwrite(Temp,fullfile(paths.seed,...
        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_epi.nii.gz',...
            MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
            num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
    end
  
% --------------------Transform images into MNI space-------------------- %

             % ----------------Whole-Brain ---------------%
% Generate names based on seed shape set by configs.
    if configs.seedShape == 0 % Sphere
        % Input Fisher's Z image in EPI space.
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_epi.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
        % Output MNI space image name.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
    elseif configs.seedShape == 1 % Cuboid
        % Input Fisher's Z image in EPI space.
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_epi.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
        % Output MNI space image name.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
    end
    
% Transformation field and matrix to MNI space.
        fileMat = fullfile(paths.EPI,'epi_2_MNI_final.mat');
        fileWarp = fullfile(paths.reg,'T12MNI_warp.nii.gz');
        
% Select reference image based on specified resolution.
        if configs.MNIres == 1
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
        elseif configs.MNIres == 2
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
        else
            error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
        end
        
% Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        
             % ------------------Seed-Only ---------------%
 % Generate names based on seed shape set by configs.
    if configs.seedShape == 0 % Sphere
        % Input seed Fisher's Z
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_epi.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
        % Output seed MNI space image.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));  
    elseif configs.seedShape == 1 % Cuboid
        % Input seed Fisher's Z
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_epi.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
        % Output seed MNI space image.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
    end
        
% Select reference image based on specified resolution.
        if configs.MNIres == 1
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
        elseif configs.MNIres == 2
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
        else
            error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
        end
        
% Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        
% --------------Copy images into group results directory----------------- %

             % ----------------Whole-Brain ---------------%
    if configs.seedShape == 0 % Sphere
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad)); 
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
    elseif configs.seedShape == 1 % Cuboid
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
    end
    sentence = sprintf('cp %s %s',fileIn,fileOut);
    [~,result]=system(sentence);
    
             % --------------- Seed-Only ----------------- %
    if configs.seedShape == 0 % Sphere
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
    elseif configs.seedShape == 1 % Cuboid
        fileIn = fullfile(paths.seed,...
            sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_MNI.nii.gz',...
        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
    end
    sentence = sprintf('cp %s %s',fileIn,fileOut);
    [~,result]=system(sentence);  

% ----------------------Smooth the MNI space data------------------------ %

                % --------convert FWHM to sigma-------- %
        configs.sigma = (configs.fwhm/sqrt(8*log(2)));
        
             % ----------------Whole-Brain ---------------%
% Generate image name. 
    if configs.seedShape == 0 % Sphere
            % Read in MNI space Fisher's Z image
        brainMNI = MRIread(fullfile(paths.resultscp,...
                sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad)));
            % Replace NaN with zeros
        brainMNI.vol(isnan(brainMNI.vol))=0;
            % Overwrite the NaN containing image with zero containing image
        MRIwrite(brainMNI,fullfile(paths.resultscp,...
                sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad)));
        % Smooth the data
        fileIn = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
        fileSmoothed = fullfile(paths.resultscp,...
            sprintf('%.0ffwhm_%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
        configs.fwhm,list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
    elseif configs.seedShape == 1 % Cuboid
        % Read in MNI space Fisher's Z image
        brainMNI = MRIread(fullfile(paths.resultscp,...
                sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
        % Replace NaN with zeros
        brainMNI.vol(isnan(brainMNI.vol))=0;
        % Overwrite the NaN containing image with zero containing image
        MRIwrite(brainMNI,fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
        % Smooth the data
        fileIn = fullfile(paths.resultscp,...
            sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
        fileSmoothed = fullfile(paths.resultscp,...
            sprintf('%.0ffwhm_%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
        configs.fwhm,list.subjects(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
    end
        
% Smooth the image.    
        sentence = sprintf('fslmaths %s -s %.8f %s',fileIn,configs.sigma,fileSmoothed);
            [~,result]=system(sentence);
    end
end
%% -----------------Use Provided Image as Seed Input-------------------- %%
if flags.runImage == 1
    if exist(seedImage,'file')
        % find dimensions of seed and reference MRI volumes        
        if configs.MNIres == 1
            refImage = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
        elseif configs.MNIres == 2
            refImage = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
        else
            error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
        end
        seedVol = MRIread(seedImage);
        refVol = MRIread(refImage);
        
        if seedVol.volsize == refVol.volsize 
            disp('ROI and Reference Volume dimensions match. Proceeding...')
        else
            sprintf('Reference volume: %s', refImage)
            sprintf('Seed volume: %s',seedImage)
            warning('Seed and Reference Volume dimensions do not match!')
            warning('Reslice Seed volume to match Reference volume. Exiting...')
            return
        end
        
        % ------------------------MNI -> T1 --------------------- %
        fileWarp = fullfile(paths.data,list.subjects(i).name,...
            'T1/registration/MNI2T1_warp.nii.gz');
        fileRef = fullfile(paths.EPI,'rT1_brain_mask.nii.gz');
        fileMat = fullfile(paths.EPI,'MNI_2_epi_final.mat');
        
        fileOut = fullfile(paths.seed,sprintf('%s_epi.nii.gz',seedLabel));
        
        % apply warp and matrix to bring seed into EPI space.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --postmat=%s --interp=nn',...
            seedImage,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        % Binarize the EPI space seed
        sentence = sprintf('fslmaths %s -bin %s',fileOut,fileOut);
        [~,result]=system(sentence);
        
        % Set seedMask as the newly generated ROI
        seedMask = MRIread(fileOut);
        if configs.GMmask == 1  % Include Only GM ROI voxels
            seedMask = logical(seedMask.vol.*subjMaskGM);
        else                    % Use ALL GM voxels
            seedMask = logical(seedMask.vol);
        end
        
        % Reshape 4D mat to a 2D, where rows=voxles and colums=volumes/time
        voxels = reshape(subjVOL,[size(subjVOL,1)*size(subjVOL,2)*...
            size(subjVOL,3),numVols]);
        % Extract time-series of seed voxels.
        seed_voxels=voxels(seedMask,:);
        % Get mean time-series of the seed
        seed=mean(seed_voxels,1);
        
        % Correlate mean seed time-series with whole-brain
        C.vector = corr(seed',voxels');
        % If single voxel native space; avoid infinite Fisher by replacing
        % corr 1 with .99 (2018.02.13 JC)
        if nnz(C.vector==1)==1
            [row, col]=find(C.vector==1);
            C.vector(row,col)=.99;
        end
        % Isolate mean seed time-series correlations with seed voxels.
        S.vector = C.vector;
        S.vector(~seedMask) = 0.;
        
        % Fisher's Z transformation of correlations
        C.vector = 0.5*log((1+C.vector)./(1-C.vector));
        S.vector = 0.5*log((1+S.vector)./(1-S.vector));
        
        % -----------Write NIfTI volumes Fisher's Z of correlations-------------- %
        
        % ----------------Whole-Brain ---------------%
        C.vol = nan(size(subjMaskGM));
        C.vol = reshape(C.vector,size(C.vol));
        Temp.vol = C.vol;
        MRIwrite(Temp,fullfile(paths.seed,...
            sprintf('conn_%s_wholebrain_epi.nii.gz',seedLabel)));
        
        % ------------------Seed-Only ---------------%
        S.vol = zeros(size(subjMaskGM));
        S.vol = reshape(S.vector,size(S.vol));
        Temp.vol = S.vol;
        MRIwrite(Temp,fullfile(paths.seed,...
            sprintf('conn_%s_seedonly_epi.nii.gz',seedLabel)));
        
        % --------------------Transform images into MNI space-------------------- %
        
        % ----------------Whole-Brain ---------------%
        % Input Fisher's Z image in EPI space.
        fileIn = fullfile(paths.seed,...
            sprintf('conn_%s_wholebrain_epi.nii.gz',seedLabel));
        % Output MNI space image name.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_%s_wholebrain_MNI.nii.gz',seedLabel));
        % Transformation field and matrix to MNI space.
        fileWarp = fullfile(paths.reg,'T12MNI_warp.nii.gz');
        fileMat = fullfile(paths.EPI,'epi_2_MNI_final.mat');
        % Select reference image based on specified resolution.
        if configs.MNIres == 1
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz';
        elseif configs.MNIres == 2
            fileRef = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
        else
            error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
        end
        % Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        
        % ------------------Seed-Only ---------------%
        % Input seed Fisher's Z
        fileIn = fullfile(paths.seed,...
            sprintf('conn_%s_seedonly_epi.nii.gz',seedLabel));
        % Output seed MNI space image.
        fileOut = fullfile(paths.seed,...
            sprintf('conn_%s_seedonly_MNI.nii.gz',seedLabel));
        % Apply transformation matrix first followed by warp field.
        sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
            fileIn,fileOut,fileRef,fileWarp,fileMat);
        [~,result]=system(sentence);
        
        % --------------Copy images into group results directory----------------- %
        
        % ----------------Whole-Brain ---------------%
        fileIn = fullfile(paths.seed,...
            sprintf('conn_%s_wholebrain_MNI.nii.gz',seedLabel));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
            list.subjects(i).name,seedLabel));
        sentence = sprintf('cp %s %s',fileIn,fileOut);
        [~,result]=system(sentence);
        
        % ------------------Seed-Only ---------------%
        fileIn = fullfile(paths.seed,...
            sprintf('conn_%s_seedonly_MNI.nii.gz',seedLabel));
        fileOut = fullfile(paths.resultscp,...
            sprintf('%s_conn_%s_seedonly_MNI.nii.gz',...
            list.subjects(i).name,seedLabel));
        sentence = sprintf('cp %s %s',fileIn,fileOut);
        [~,result]=system(sentence);
        
        % ----------------------Smooth the MNI space data------------------------ %
        
        % --------convert FWHM to sigma-------- %
        configs.sigma = (configs.fwhm/sqrt(8*log(2)));
        
        % ----------------Whole-Brain ---------------%
        % Read in MNI space Fisher's Z image
        brainMNI = MRIread(fullfile(paths.resultscp,...
            sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
            list.subjects(i).name,seedLabel)));
        % Replace NaN with zeros
        brainMNI.vol(isnan(brainMNI.vol))=0;
        % Overwrite the NaN containing image with zero containing image
        MRIwrite(brainMNI,fullfile(paths.resultscp,...
            sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
            list.subjects(i).name,seedLabel)));
        % Smooth the data
        fileIn = fullfile(paths.resultscp,...
            sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
            list.subjects(i).name,seedLabel));
        fileSmoothed = fullfile(paths.resultscp,...
            sprintf('%.0ffwhm_%s_conn_%s_wholebrain_MNI.nii.gz',...
            configs.fwhm,list.subjects(i).name,seedLabel));
        % Smooth the Images
        sentence = sprintf('fslmaths %s -s %.8f %s',fileIn,configs.sigma,fileSmoothed);
        [~,result]=system(sentence);
    else
        error('Seed image file does not exist or is not specified. Terminating...')
    end
end
end










