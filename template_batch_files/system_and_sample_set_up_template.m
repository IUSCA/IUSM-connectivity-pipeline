%                       SYSTEM_AND_SAMPLE_SET_UP
%    This contains all general paths and variables necessary to run the
%  pipeline batch files. You may edit this as necessary depending on your
%            program paths and the subjects you want to run. 
%
% Contributors:
%        Evgeny Chumin, Indiana University School of Medicine, 2018
%                       Indiana University Bloomington, 2020                    
%        John West, Indiana University School of Medicine, 2018
%%
            %------------------------------------------------%
            %  SET PATH TO THE CONNECTOME SCRIPTS DIRECTORY  %
            %------------------------------------------------%

    % Add path to connectome scripts directory
paths.scripts = '/N/dc2/projects/connectivitypipeline/example_for_andrea/jenya/IUSM-connectivity-pipeline';
addpath(genpath(paths.scripts));

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
    % ants brain extraction path
[~,antspath]=system('which antsBrainExtraction.sh');
paths.ANTS = extractBefore(antspath,'/antsBrainExtraction.sh');
clear antspath
    % Path to MRtrix
[~,mrtrix]=system('which dwi2response');
paths.MRtrix = extractBefore(mrtrix,'/dwi2response');
    % path to python (for ICA-AROMA EPI clean-up)
    % 
    % consult readme and manual in the ICA-AROMA directory in
    % connectome_scripts for instalation instructions.
%[~,paths.python]=system('which python');
    % if you want to use a different python distribution (other than that
    % on your system, commento ut line 40, and uncomment and set path in
    % line 43.
 paths.python = '/N/u/echumin/Carbonate/miniconda3/bin/python3';
%%
                    %------------------------------%
                    %  SELECT SUBJECT DIRECTORIES  %
                    %------------------------------%
    % Set the path to the directory containing you subjects.
paths.data = '/N/project/kbase-imaging/all_prelim_data';

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
%subjectList(1).name = 'SU0360'; 
%subjectList(end+1).name = 'SUBEJCT2'; % copy this line for additional subjects

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

% ALTERNTIVE if multiple spin echo fielmap scans are available
% SUBJECT1 -- T1 -- DICOMS
%          |
%          -- EPI(#) -- DICOMS (May have multiple EPI scans)
%          |         |
%          |         |      (SPIN-ECHO)     
%          -- UNWARP1 -- SEFM_AP_DICOMS 
%          |                   | 
%          |                   -- SEFM_PA_DICOMS (OR) GREFM_PHASE_DICOMS
%          -- UNWARP2 -- SEFM_AP_DICOMS 
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
configs.name.grefmFolder = 'GREFM_GUST'; % Reserved for Field Mapping series
    configs.name.GREmagdcm = 'MAG_DICOMS'; % Gradient echo FM magnitude series
    configs.name.GREphasedcm = 'PHASE_DICOMS'; % Gradient echo FM phase map series
configs.name.DWI = 'DWI';
    configs.name.unwarpFolder = 'UNWARP';
        configs.name.dcmPA = 'B0_PA_DCM'; %b0 opposite phase encoding

configs.name.dcmFolder = 'DICOMS';

%%
