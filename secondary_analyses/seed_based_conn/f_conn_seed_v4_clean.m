clearvars

%% --------------------- USER DEFINED INPUT DATA ----------------------- %%

% Set path to connectivity pipeline directory
paths.scripts=...
    '/N/u/echumin/Carbonate/IUSM-connectivity-pipeline';

% Data directory (contains subject directories)
paths.data=...
    '/N/dc2/projects/brainconnectomics/IADC-IMAS-image-processing/datadir';

% Define a path and directory name to put the results
%   If one does not exist, it will be created, if it does exist, be aware
%   that files may be overwritten.
paths.results =...
    '/N/dc2/projects/brainconnectomics/IADC-IMAS-image-processing/test_out';
%-------------------------------------------------------------------------%
% Provide the names of subdirectories used for processing
    % fMRI subdirectory name 
        epiDir = 'EPI';                 % e.g. EPI, REST, EPI1, etc.
    % Nuisance regression strategy used
        nuiDir = 'AROMA/aCompCorr';     % Some possible options:
                                        % AROMA/aCompCorr
                                        % HMPreg/PhysReg
                                        % refer to your data for subdirectories
    % String of nuisance parameter selections.
        % Refer to the TimeSeries* subdirectory for your subjects.
        nui_str = 'aroma_Gs4_pca5';
%-------------------------------------------------------------------------%
% Set up a subjectt list (several options) uncomment and edit the desired option.
    % option 1: use all subject directorties in data directories
         subjectList = [];
        
    % option 2: provide a subject list as a matlab structure named subjectList, 
    % with a .name field, saved in a .mat file (input both path and filename)
        % subjectList =...
        %    '/N/dc2/projects/brainconnectomics/IADC-IMAS-image-processing/subjectLists/subjectList_test.mat';
%-------------------------------------------------------------------------%
% if EPI-to-MNI transformation matrices already exist (previous seed),
% should they be recalculated for subsequent seeds?
% 1 = Yes (slows down the code but helpful if starting from scratch)
% 0 = No (calculate once, then use those matrices for subsequent seeds)
configs.TransMatOverwrite = 1;

%% --------------------------- DEFINE SEEDS ---------------------------- %%
% There are several options for defining seeds (1) from parcellations used
% in the pipeline, (2) MNI coortidate list and size for spherical or
% rectangular ROI placement, and (3) an MNI space binary image.
% WHICH SEED OPTIONS ARE USED IS SET IN THE CONFIGS SECTION BELOW.
%-------------------------------------------------------------------------%
% Option1:
% Provide parcellation indices;
    list.regions = [301];
    parc_name = 'schaefer300_yeo17'; % match the parcs.plabel.name set in 
                                     % the pipeline batch file
    configs.erodeLevel = 2; % this is the level of erosion on seeds. 
                            % 1 does nothing.
                            % 3 most likely recommended. 5 may remove seeds entirely.
                            % even numbers will produce uneven erosions, but are possible.
%-------------------------------------------------------------------------%
% Option2:
% Provide MNI coordinates X Y Z (tab delimited) 3-column test file as origins
% for generation of seeds.
    seedCoordList = ...
        '/N/u/echumin/Carbonate/r_caud_seed.txt';

% Define Seed parameters:   
    configs.seedShape =0;   % =0 Sphere
                            % =1 Rectangle
% Define Seed Size
    configs.seedRad = 5;        % Radius in mm of the sphere
    configs.seedBox = [5 5 5];  % Must be odd values!
                                % [ X   Y   Z] dimentions IN VOXELS of the box.
                                % [L-R A-P S-I]
%-------------------------------------------------------------------------%
% Option3:
% Provide a path to an image you wish to use as a seed, as well as a label
% for the seed.
% This image will binarized in the code, therefore if you wish to
% impose a threshold, do so before providing image as input.
% The input must be in MNI152_T1 space.
seedLabel = 'L_ins';
seedImage = ...
    '/N/u/echumin/Carbonate/R_insula_MNI_1mm.nii.gz';

%% -------------------------- OTHER CONFIGS ---------------------------- %%

% Use only Gray Matter voxels 
% (applicable only to seed coordinate and image methods)
configs.GMmask = 1; % 1 = "On"; 0 = "Off"

% Number of volumes dropped off the start and end of EPI series.
configs.dropTR = 5; % defaults to 15

% Specify the resolution of MNI space seedconn results: 1mm or 2mm.
% For image seed input this should also match the resolution of seed image.
configs.MNIres = 1;

% Smoothing kernel Full-width Half-max (mm)
configs.fwhm = 4;

%% ------------------- CHOOSE WHICH INPUT TO USE ----------------------- %%
% Shen region IDs as input
flags.runParc = 1;

% List of MNI coordinates as input
flags.runMNIcoords = 1;

% MNI space image as input
flags.runImage = 1; % Currently being tested

%% ----------------------------------------------------------------------%%

% Latest revistion: May 2020, Jenya Chumin

%% --------------ALL CODE BELOW THIS LINE IS NO TOUCHIE!---------------- %%
% Finding FSL path
[~,fpath]=system('which fsl'); fpath = extractBefore(fpath,'/bin');
paths.FSL = [fpath '/bin'];

% Addings connectome spripts path (Includes NIfTI toolbox)
addpath(genpath(paths.scripts))

if ~exist(paths.results,'dir')
    mkdir(paths.results)
else
    fprintf(2,'Results will be written to existing directory:\n %s\n',paths.results)
end

% Generate a list of subjects based on provided input.
if isempty(subjectList)
    disp('SubjectList variable empty. Using all directories in data path.')
    subjectList=dir(paths.data);
    subjectList(1:2)=[];
else
    load(subjectList)
    if isstruct(subjectList)
        disp('Using loaded subjectList')
    else
        disp('No subjectList variable in file OR it is not a structure.')
        disp('Exiting....')
        return
    end
end
%-------------------------------------------------------------------------%
paths.allOUT = fullfile(paths.results,strcat(epiDir,'_',nui_str));
    if ~exist(paths.allOUT,'dir')
        mkdir(paths.allOUT)
    end
 numSubj=length(subjectList);
% ---------------------Looping Across Subjects--------------------------- %
for i=1:numSubj     
    paths.EPI=fullfile(paths.data,subjectList(i).name,epiDir);
    if exist(paths.EPI,'dir')
        fprintf('Running: %d/%d %s\n',i,numSubj,subjectList(i).name)
        % Set subject specific paths
        paths.GSPCA = fullfile(paths.EPI,nuiDir);
        paths.seed = fullfile(paths.GSPCA,'SeedConnectivity');
        if ~exist(paths.seed,'dir')
            mkdir(paths.seed)
        else
            fprintf(2,'Results for %s will be written to existing directory:\n %s\n',...
                subjectList(i).name,paths.seed)
        end
        % Read in the EPI 4D volume for the subject
        fileNii = fullfile(paths.GSPCA,strcat('7_epi_',nui_str,'.nii.gz'));
        if exist(fileNii,'file')
            subjNii = MRIread(fileNii);
        else
            fprintf(2,'Processed epi data file does not exist.\n')
            fprintf(2,'Terminating subject: %s\n',subjectList(i).name)
            return
        end
        % Drop early and late volumes
        subjVOL = subjNii.vol(:,:,:,configs.dropTR+1:end-configs.dropTR);clear subjNii;
        % Calculate numver of remaining time points
        numVols = size(subjVOL,4);
            % add if no parc_name given use GM mask; if no GM mask quit
        % Read in the EPI space parcellation
        if isempty(parc_name) || ~exist('parc_name','var')
            fileGM=fullfile(paths.EPI,'rT1_GM_mask.nii.gz');
            if exist(fileGM,'file')
                subjGM = MRIread(fileGM);
                subjMaskGM = subjGM.vol > 0;
                Temp = subjGM;
            else
                fprintf(2,'Parcellation name not specified and\n')
                fprintf(2,'no epi GM mask found. Terminating subject...\n')
                return
            end
        else
            filePARC=fullfile(paths.EPI,strcat('rT1_GM_parc_',parc_name,'_clean.nii.gz'));
            if exist(filePARC,'file')
                subjPARC = MRIread(filePARC);
                % Create mask from parcellation
                subjMaskGM = subjPARC.vol > 0;
                % Create a Temp volume structure for writing out NIfTI files
                Temp = subjPARC;
            else
                fprintf(2,'parc_name provided does not match a parcellation image.\n')
                fileGM=fullfile(paths.EPI,'rT1_GM_mask.nii.gz');
                if exist(fileGM,'file')
                    subjGM = MRIread(fileGM);
                    subjMaskGM = subjGM.vol > 0;
                    Temp = subjGM;
                else
                    fprintf(2,'no epi GM mask found. Terminating subject...\n')
                    return
                end
            end
        end
% -----------Prepare transformation matrices for MNI space--------------- %

        % EPI -> T1 linear dof6 and bbr (inverse of T1 -> EPI)
        fileMat = fullfile(paths.EPI,'T1_2_epi_dof6_bbr.mat');
        fileMatInv = fullfile(paths.EPI,'epi_2_T1_bbr_dof6.mat');
        if ~exist(fileMatInv,'file')
            if exist(fileMat,'file')
                disp('Generating an epi 2 T1 transformation')
                sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',paths.FSL,fileMatInv,fileMat);
                [status,result] = system(sentence);
            else
                fprintf(2,'Transformations between T1 and EPI not found.\n')
                fprintf(2,'Terminating subject: %s\n',subjectList(i).name)
                return
            end
        else
            fprintf('Will use existing EPI to T1 transformation.\n')
        end

        % Combine with T12MNI_dof6
        paths.reg = fullfile(paths.data,subjectList(i).name,'T1/registration');
        fileMat1 = fullfile(paths.reg,'T12MNI_dof6.mat');
        fileMatJoint = fullfile(paths.EPI,'epi_2_MNI_dof6.mat');
        if ~exist(fileMatJoint,'file') || configs.TransMatOverwrite == 1
            sentence = sprintf('%s/convert_xfm -omat %s -concat %s %s',paths.FSL,fileMatJoint,fileMat1,fileMatInv);
            [status,result] = system(sentence);
        else
            fprintf('Using available epi_2_MNI_dof6\n')
        end

        % Finally, apply T12MNI_dof12
        fileMat1 = fullfile(paths.reg,'T12MNI_dof12.mat');
        fileMat2 = fullfile(paths.EPI,'epi_2_MNI_dof6.mat');
        fileMatJoint2 = fullfile(paths.EPI,'epi_2_MNI_final.mat');
        if ~exist(fileMatJoint2,'file') || configs.TransMatOverwrite == 1
            sentence = sprintf('%s/convert_xfm -omat %s -concat %s %s',paths.FSL,fileMatJoint2,fileMatJoint,fileMat2);
            [status,result] = system(sentence);
        else
            fprintf('Using available epi_2_MNI_final\n')
        end
        % Make an inverse
        fileMatInv2 = fullfile(paths.EPI,'MNI_2_epi_final.mat');
        if ~exist(fileMatInv2,'file') || configs.TransMatOverwrite == 1
            sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',paths.FSL,fileMatInv2,fileMatJoint2);
            [status,result] = system(sentence);
        else
            fprintf('Using available MNI_2_epi_final\n')
        end
   % T1 -> MNI nonlinear warp  - generated in the pipeline.
%-------------------------------------------------------------------------%

%% ------------------Looping Across Shen Seeds-------------------------- %%
        if flags.runParc == 1
            for j=1:length(list.regions)
                fprintf('Running %s region number %d\n',parc_name,list.regions(j))
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
                C.vector(C.vector>=1) = 0.99999; % make it well-behaved (C <= 1)
                % Isolate mean seed time-series correlations with seed voxels.
                S.vector = C.vector;
                S.vector(~seedMask) = 0.;

                % Fisher's Z transformation of correlations
                C.vector = 0.5*log((1+C.vector)./(1-C.vector));
                S.vector = 0.5*log((1+S.vector)./(1-S.vector));

% -----------Wrie NIfTI volumes Fisher's Z of correlations-------------- %

                % ----------------Whole-Brain ---------------%
                C.vol = nan(size(subjMaskGM));
                C.vol = reshape(C.vector,size(C.vol));
                Temp.vol = C.vol;
                MRIwrite(Temp,fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_wholebrain_epi.nii.gz',parc_name,list.regions(j))));

                % ------------------Seed-Only ---------------%
                S.vol = zeros(size(subjMaskGM));
                S.vol = reshape(S.vector,size(S.vol));
                Temp.vol = S.vol;
                MRIwrite(Temp,fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_seedonly_epi.nii.gz',parc_name,list.regions(j))));

% --------------------Transform images into MNI space-------------------- %

                % ----------------Whole-Brain ---------------%
                % Input Fisher's Z image in EPI space.
                fileIn = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_wholebrain_epi.nii.gz',parc_name,list.regions(j)));
                % Output MNI space image name.
                fileOut = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',parc_name,list.regions(j)));
                % Transformation field and matrix to MNI space.
                fileMat = fullfile(paths.EPI,'epi_2_MNI_final.mat');
                fileWarp = fullfile(paths.reg,'T12MNI_warp.nii.gz');
                % Select reference image based on specified resolution.
                if configs.MNIres == 1
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_1mm.nii.gz');
                elseif configs.MNIres == 2
                    fileRef = fullfile(fmath,'data/standard/MNI152_T1_2mm.nii.gz');
                else
                    error('Error. \nInvalid resolution of MNI template: %.0f',configs.MNIres)
                end
                % Apply transformation matrix first followed by warp field.
                sentence = sprintf('%s/applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
                    paths.FSL,fileIn,fileOut,fileRef,fileWarp,fileMat);
                [~,result]=system(sentence);

                % ------------------Seed-Only ---------------%
                % Input seed Fisher's Z
                fileIn = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_seedonly_epi.nii.gz',parc_name,list.regions(j)));
                % Output seed MNI space image.
                fileOut = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_seedonly_MNI.nii.gz',parc_name,list.regions(j)));
                % Apply transformation matrix first followed by warp field.
                sentence = sprintf('applywarp -i %s -o %s -r %s -w %s --premat=%s --interp=nn',...
                    fileIn,fileOut,fileRef,fileWarp,fileMat);
                [~,result]=system(sentence);

                % --------------Copy images into group results directory----------------- %

                % ----------------Whole-Brain ---------------%
                fileIn = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',parc_name,list.regions(j)));
                fileOut = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,parc_name,list.regions(j)));
                sentence = sprintf('cp %s %s',fileIn,fileOut);
                [~,result]=system(sentence);

                % ------------------Seed-Only ---------------%
                fileIn = fullfile(paths.seed,...
                    sprintf('conn_%s_seed%03d_eroded_seedonly_MNI.nii.gz',parc_name,list.regions(j)));
                fileOut = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seed%03d_eroded_seedonly_MNI.nii.gz',...
                    subjectList(i).name,parc_name,list.regions(j)));
                sentence = sprintf('cp %s %s',fileIn,fileOut);
                [~,result]=system(sentence);

                % ----------------------Smooth the MNI space data------------------------ %

                % --------convert FWHM to sigma-------- %
                configs.sigma = (configs.fwhm/sqrt(8*log(2)));

                % ----------------Whole-Brain ---------------%
                % Read in MNI space Fisher's Z image
                brainMNI = MRIread(fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,parc_name,list.regions(j))));
                %         % Replace NaN with zeros
                brainMNI.vol(isnan(brainMNI.vol))=0;
                %         % Overwrite the NaN containing image with zero containing image
                MRIwrite(brainMNI,fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,parc_name,list.regions(j))));
                %         % Smooth the data
                fileIn = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,parc_name,list.regions(j)));
                fileSmoothed = fullfile(paths.allOUT,...
                    sprintf('%.0ffwhm_%s_conn_%s_seed%03d_eroded_wholebrain_MNI.nii.gz',...
                    configs.fwhm,subjectList(i).name,parc_name,list.regions(j)));
                sentence = sprintf('%s/fslmaths %s -s %.8f %s',paths.FSL,fileIn,configs.sigma,fileSmoothed);
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
            fileSource = fullfile(fpath,'data/standard/MNI152_T1_1mm_brain.nii.gz'); % used by the pipeline - DO NOT CHANGE
            fileDest = fullfile(paths.data,subjectList(i).name,...
                'T1/registration/T1_dof12.nii.gz');
            filewarp = fullfile(paths.data,subjectList(i).name,...
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
            fileSource = fullfile(paths.data,subjectList(i).name,...
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
                sentence = sprintf('%s/fslmaths %s -roi %s 1 %s 1 %s 1 0 1 %s'...
                    ,paths.FSL,ImageIN...
                    ,num2str(voxelCoords(l,1)),num2str(voxelCoords(l,2)),num2str(voxelCoords(l,3))...
                    ,SeedOut);
                [~,result]=system(sentence);

                % ---------------------Growing the seeds-------------------- %
                if configs.seedShape == 0 % Sphere
                    SeedName2 = sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphereROI.nii.gz',...
                        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3));
                    SeedOut2 = fullfile(paths.seed,SeedName2);
                    sentence = sprintf('%s/fslmaths %s -kernel sphere %d -fmean -bin %s',...
                        paths.FSL,SeedOut,configs.seedRad,SeedOut2);
                    [~,result]=system(sentence);

                else % Rectangle
                    SeedName2 = sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboidROI.nii.gz',...
                        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3));
                    SeedOut2 = fullfile(paths.seed,SeedName2);
                    sentence = sprintf('%s/fslmaths %s -kernel boxv3 %s %s %s -fmean -bin %s',...
                        paths.FSL,SeedOut,...
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
                C.vector(C.vector>=1) = 0.99999; % make it well-behaved (C <= 1)
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
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_1mm.nii.gz');
                elseif configs.MNIres == 2
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_2mm.nii.gz');
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
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_1mm.nii.gz');
                elseif configs.MNIres == 2
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_2mm.nii.gz');
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
                    fileOut = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
                elseif configs.seedShape == 1 % Cuboid
                    fileIn = fullfile(paths.seed,...
                        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
                    fileOut = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
                end
                sentence = sprintf('cp %s %s',fileIn,fileOut);
                [~,result]=system(sentence);

                % --------------- Seed-Only ----------------- %
                if configs.seedShape == 0 % Sphere
                    fileIn = fullfile(paths.seed,...
                        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_MNI.nii.gz',...
                        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
                    fileOut = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_seedonly_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
                elseif configs.seedShape == 1 % Cuboid
                    fileIn = fullfile(paths.seed,...
                        sprintf('conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_MNI.nii.gz',...
                        MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
                    fileOut = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_seedonly_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
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
                    brainMNI = MRIread(fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad)));
                    % Replace NaN with zeros
                    brainMNI.vol(isnan(brainMNI.vol))=0;
                    % Overwrite the NaN containing image with zero containing image
                    MRIwrite(brainMNI,fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad)));
                    % Smooth the data
                    fileIn = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
                    fileSmoothed = fullfile(paths.allOUT,...
                        sprintf('%.0ffwhm_%s_conn_seed_MNI_%.0f_%.0f_%.0f_sphere_%.0fmmRad_wholebrain_MNI.nii.gz',...
                        configs.fwhm,subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),configs.seedRad));
                elseif configs.seedShape == 1 % Cuboid
                    % Read in MNI space Fisher's Z image
                    brainMNI = MRIread(fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
                    % Replace NaN with zeros
                    brainMNI.vol(isnan(brainMNI.vol))=0;
                    % Overwrite the NaN containing image with zero containing image
                    MRIwrite(brainMNI,fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3)))));
                    % Smooth the data
                    fileIn = fullfile(paths.allOUT,...
                        sprintf('%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
                    fileSmoothed = fullfile(paths.allOUT,...
                        sprintf('%.0ffwhm_%s_conn_seed_MNI_%.0f_%.0f_%.0f_cuboid_%s_%s_%sDims_wholebrain_MNI.nii.gz',...
                        configs.fwhm,subjectList(i).name,MNIcoords(l,1),MNIcoords(l,2),MNIcoords(l,3),...
                        num2str(configs.seedBox(1,1)),num2str(configs.seedBox(1,2)),num2str(configs.seedBox(1,3))));
                end

                % Smooth the image.
                sentence = sprintf('%s/fslmaths %s -s %.8f %s',paths.FSL,fileIn,configs.sigma,fileSmoothed);
                [~,result]=system(sentence);
            end
        end
        %% -----------------Use Provided Image as Seed Input-------------------- %%
        if flags.runImage == 1
            if exist(seedImage,'file')
                % find dimensions of seed and reference MRI volumes
                if configs.MNIres == 1
                    refImage = fullfile(fpath,'data/standard/MNI152_T1_1mm.nii.gz');
                elseif configs.MNIres == 2
                    refImage = fullfile(fpath,'data/standard/MNI152_T1_2mm.nii.gz');
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
                fileWarp = fullfile(paths.data,subjectList(i).name,...
                    'T1/registration/MNI2T1_warp.nii.gz');
                fileRef = fullfile(paths.EPI,'rT1_brain_mask.nii.gz');
                fileMat = fullfile(paths.EPI,'MNI_2_epi_final.mat');

                fileOut = fullfile(paths.seed,sprintf('%s_epi.nii.gz',seedLabel));

                % apply warp and matrix to bring seed into EPI space.
                sentence = sprintf('%s/applywarp -i %s -o %s -r %s -w %s --postmat=%s --interp=nn',...
                    paths.FSL,seedImage,fileOut,fileRef,fileWarp,fileMat);
                [~,result]=system(sentence);
                % Binarize the EPI space seed
                sentence = sprintf('%s/fslmaths %s -bin %s',paths.FSL,fileOut,fileOut);
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
                C.vector(C.vector>=1) = 0.99999; % make it well-behaved (C <= 1)
                % If single voxel native space; avoid infinite Fisher by replacing
                % corr 1 with .99 (2018.02.13 JC)
    %                 if nnz(C.vector==1)==1
    %                     [row, col]=find(C.vector==1);
    %                     C.vector(row,col)=.99;
    %                 end
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
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_1mm.nii.gz');
                elseif configs.MNIres == 2
                    fileRef = fullfile(fpath,'data/standard/MNI152_T1_2mm.nii.gz');
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
                fileOut = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,seedLabel));
                sentence = sprintf('cp %s %s',fileIn,fileOut);
                [~,result]=system(sentence);

                % ------------------Seed-Only ---------------%
                fileIn = fullfile(paths.seed,...
                    sprintf('conn_%s_seedonly_MNI.nii.gz',seedLabel));
                fileOut = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_seedonly_MNI.nii.gz',...
                    subjectList(i).name,seedLabel));
                sentence = sprintf('cp %s %s',fileIn,fileOut);
                [~,result]=system(sentence);

                % ----------------------Smooth the MNI space data------------------------ %

                % --------convert FWHM to sigma-------- %
                configs.sigma = (configs.fwhm/sqrt(8*log(2)));

                % ----------------Whole-Brain ---------------%
                % Read in MNI space Fisher's Z image
                brainMNI = MRIread(fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,seedLabel)));
                % Replace NaN with zeros
                brainMNI.vol(isnan(brainMNI.vol))=0;
                % Overwrite the NaN containing image with zero containing image
                MRIwrite(brainMNI,fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,seedLabel)));
                % Smooth the data
                fileIn = fullfile(paths.allOUT,...
                    sprintf('%s_conn_%s_wholebrain_MNI.nii.gz',...
                    subjectList(i).name,seedLabel));
                fileSmoothed = fullfile(paths.allOUT,...
                    sprintf('%.0ffwhm_%s_conn_%s_wholebrain_MNI.nii.gz',...
                    configs.fwhm,subjectList(i).name,seedLabel));
                % Smooth the Images
                sentence = sprintf('fslmaths %s -s %.8f %s',fileIn,configs.sigma,fileSmoothed);
                [~,result]=system(sentence);
            else
                error('Seed image file does not exist or is not specified. Terminating...')
            end
        end
    else
        fprintf(2,'No %s directory found for %s. Skipping subject...\n',epiDir,subjectList(i).name)
    end
end










