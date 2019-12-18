%                       SYSTEM_AND_SAMPLE_SET_UP
%    This contains all general paths and variables necessary to run the
%  pipeline batch files. You may edit this as necessary depending on your
%            program paths and the subjects you want to run. 
%
% Contributors:
%        Evgeny Chumin, Indiana University School of Medicine, 2018
%        John West, Indiana University School of Medicine, 2018
%%
            %------------------------------------------------%
            %  SET PATH TO THE CONNECTOME SCRIPTS DIRECTORY  %
            %------------------------------------------------%

    % Add path to connectome scripts directory
paths.scripts = '/home/echumin/Documents/PROJECTS/image_processing/IUSM-connectivity-pipeline/connectome_scripts';
addpath(paths.scripts);
    % path to use MRIread MRIwrite
addpath(fullfile(paths.scripts,'toolbox_matlab_nifti'));
    % path to templates in MNI
paths.MNIparcs =fullfile(paths.scripts,'templates/MNIparcs');
    % path to T1 denoiser
addpath(genpath(fullfile(paths.scripts,'MRIDenoisingPackage')));

%%  (This may/should already be set in your .bashrc)
    % path to FSL bin directory
[~,fpath]=system('which fsl'); fpath = extractBefore(fpath,'/bin');
paths.FSL = [fpath '/bin'];
    % FSL setup
FSLsetup = ['FSLDIR=' fpath '; . ${FSLDIR}/etc/fslconf/fsl.sh; PATH=${FSLDIR}/bin:${PATH}; export FSLDIR PATH'];
clear fpath
    % Path to AFNI
[~,apath]=system('which afni');
paths.AFNI = apath(1:end-6);
clear apath
    % Path to MRIcroGL
paths.MRIcroGL = '/N/soft/rhel6/mricrogl/20170808/mricrogl_lx';
    % ants brain extraction path
[~,antspath]=system('which ants.sh');
paths.ANTS = extractBefore(antspath,'/ants.sh');
clear antspath
    % path to python (for ICA-AROMA EPI clean-up)
    % currently using python2.7 in the miniconda environment
    % consult readme and manual in the ICA-AROMA directory in
    % connectome_scripts for instalation instructions.
paths.python = '~/anaconda3/bin/python3.7';
%%
                    %------------------------------%
                    %  SELECT SUBJECT DIRECTORIES  %
                    %------------------------------%
    % Set the path to the directory containing you subjects.
paths.data = '/home/echumin/Documents/PROJECTS/image_processing/datadir/';

% NOTE: For supercomputing job submissions DO NOT specify a subjectList
% here. It is generated separately by the PBS job generator. 
                    if ~exist('subjectList','var')
                        %{
                        If a subjectList variable is 
                        already in the matlab environment
                        the code below will not execute.
                        This prevents conflict with 
                        supercomputing submissions.
                        %}
    % generate a list of subjects from directories in path
subjectList =dir(paths.data); subjectList(1:2)=[]; %#ok<*NASGU> %remove '.' and '..'

    % If you wish to exclude subjects from the above generated list, use
    % the below line, replacing SUBJECT1 with the subject you want to
    % exclude. Copy and paste the line several times to exclude multiple
    % subjects.
    
% idx=find(strcmp({subjectList.name},'SUBJECT1')==1); subjectList(idx)=[];

    % If you only wish to process a specific subject or set of subjects,
    % use the following three lines as example. If processing more that 2
    % subjects copy and paste the second line as necessary.
    
%clear subjectList %remove the above generated list
%subjectList(1).name = '10004_6'; % copy this line for additional subjects

                    end

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
%          |                   | 
%          |                   -- SEFM_PA_DICOMS (OR) GREFM_PHASE_DICOMS
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
configs.name.DWI = 'DWI';
    configs.name.unwarpFolder = 'UNWARP';
        configs.name.dcmPA = 'B0_PA_DCM'; %b0 opposite phase encoding

configs.name.dcmFolder = 'DICOMS';
configs.name.dcmFiles = 'IMA'; % Dicom file extension
configs.name.niiFiles = 'nii'; % Nifti-1 file extension

%%
