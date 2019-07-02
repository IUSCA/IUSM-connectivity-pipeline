%                       SYSTEM_AND_SAMPLE_SET_UP
%    This contains all general paths and variables necessary to run the
%  pipeline batch files. You may edit this as necessary depending on your
%            program paths and the subjects you want to run. 
%
% Contributors:
%        Evgeny Chumin, Indiana University School of Medicine, 2018
%%
            %------------------------------------------------%
            %  SET PATH TO THE CONNECTOME SCRIPTS DIRECTORY  %
            %------------------------------------------------%
    % Add path to connectome scripts directory
paths.scripts = '/usr/local/IUSM-connectivity-pipeline/connectome_scripts';
addpath(paths.scripts);
    % path to use MRIread MRIwrite
addpath(fullfile(paths.scripts,'toolbox_matlab_nifti/'));
    % path to templates in MNI
paths.MNIparcs =fullfile(paths.scripts,'templates/MNIparcs');
    % path to T1 denoiser
addpath(genpath(fullfile(paths.scripts,'/MRIDenoisingPackage')));

%%  (This may/should already be set in your .bashrc)
    % path to FSL bin directory
paths.FSL = '/usr/local/fsl/bin';
    % FSL setup
FSLsetup = 'FSLDIR=/usr/local/fsl; . ${FSLDIR}/etc/fslconf/fsl.sh; PATH=${FSLDIR}/bin:${PATH}; export FSLDIR PATH';
%FSLsetup = 'FSLDIR=/data04/Zikai/IUSM-connectivity-pipeline/fsl/; . ${FSLDIR}/etc/fslconf/fsl.sh; PATH=${FSLDIR}/bin:${PATH}; export FSLDIR PATH';

    % Path to feat
paths.feat = '/usr/local/fsl/bin/feat';
    % Path to AFNI
paths.AFNI = '/usr/local/afni';
    % Path to MRIcroGL
paths.MRIcroGL = '/usr/local/mricrogl';
    % Camino setup
paths.CaminoSetup=sprintf('PATH=%s:${PATH}',fullfile('/usr/local/camino','bin'));
    % CaminoTrackVis setup
paths.CamTrackSetup=sprintf('PATH=%s:${PATH}',fullfile('/usr/local/camino-trackvis','bin'));
    % DTItk setup
paths.DTItkSetup=sprintf('PATH=%s:${PATH}',fullfile('/usr/local/dtitk','bin'));

%% ICA-AROMA paths set up
    % path to ICA-AROMA 
paths.aroma = '/usr/local/fsl/ICA-AROMA/ICA_AROMA.py'; % the program of ica-aroma has to be a python files
    % path to standard images 
paths.stdImg = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain'; % paths to standard image
    % ICA-AROMA directory name (optional, if ICA-AROMA has been processed
    % prior to running to pipeline
configs.name.ica_aroma_folder = 'ICA_AROMA';

%%
                    %------------------------------%
                    %  SELECT SUBJECT DIRECTORIES  %
                    %------------------------------%
    % Set the path to the directory containing you subjects.
paths.data = '/XXXX/CONNECTIVITY/datadir/';
    % generate a list of subjects from directories in path
%subjectList =dir(paths.data); subjectList(1:2)=[]; %#ok<*NASGU> %remove '.' and '..'

    % If you wish to exclude subjects from the above generated list, use
    % the below line, replacing SUBJECT1 with the subject you want to
    % exclude. Copy and paste the line several times to exclude multiple
    % subjects.
    
% idx=find(strcmp({subjectList.name},'SUBJECT1')==1); subjectList(idx)=[];

    % If you only wish to process a specific subject or set of subjects,
    % use the following three lines as example. If processing more that 2
    % subjects copy and paste the second line as necessary.
    
clear subjectList %remove the above generated list
subjects = ["NNNN0001"; "NNNN0002"];

% A more convenient  way for user to define subjectList
for i = 1:length(subjects)
   subjectList(i).name = char(subjects(i)); 
end

%subjectList(end+1).name = 'SUBJECT2'; % copy this line for additional subjects

%%
                    %------------------------------%
                    %  SET UP DIRECTORY STRUCTURE  %
                    %------------------------------%
% The following diagrapm is a sample directory tree for a single subject.
% Following that are configs you can use to set your own names if different
% from sample structure.

% SUBJECT1 -- T1 -- DICOMS
%          |
%          -- EPI(#) -- DICOMS (May have multiple EPI scans)
%          |         |
%          |         |               (SPIN-ECHO)       (GRADIENT ECHO)
%          |         -- UNWARP -- SEFM_AP_DICOMS (OR) GREFM_MAG_DICOMS
%          |         |         | 
%          |         |         -- SEFM_PA_DICOMS (OR) GREFM_PHASE_DICOMS
%          |         |
%          |         -- UNWARPED uf*.nii.gz (MELODIC UNWARPED IMAGES)
%          |         |
%          |         |
%          |         -- FEAT -- FEAT_PREP
%          |
%          -- DWI -- DICOMS
%                 |
%                 -- UNWARP -- B0_PA_DCM

configs.name.T1 = 'T1';
configs.name.epiFolder = 'EPI';
    configs.name.sefmFolder = 'UNWARP'; % Reserved for Field Mapping series
        configs.name.APdcm = 'SEFM_AP_DICOMS'; % Spin Echo A-P
        configs.name.PAdcm = 'SEFM_PA_DICOMS'; % Spin Echo P-A
        configs.name.GREmagdcm = 'GREFM_MAG_DICOMS'; % Gradient echo FM magnitude series
        configs.name.GREphasedcm = 'GREFM_PHASE_DICOMS'; % Gradient echo FM phase map series
    configs.name.melodicUnwarpedFolder = 'UNWARPED'; % Unwarped images (from Melodic)
configs.name.DWI = 'DWI';
    configs.name.unwarpFolder = 'UNWARP';
        configs.name.dcmPA = 'B0_PA_DCM'; % b0 opposite phase encoding

configs.name.dcmFolder = 'DICOMS';
configs.name.dcmFiles = 'dcm'; % Dicom file extension
configs.name.niiFiles = 'nii'; % Nifti-1 file extension

%%