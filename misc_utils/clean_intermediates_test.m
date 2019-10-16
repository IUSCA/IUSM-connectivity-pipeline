%               INTERMEDIATE FILE COMPRESSION/CLEANUP
%                       clean_intermediates_test
%  This code is mean to be used after all the processing and QA has been
%  done. It will compress the majority of .nii and .nii.gz into a single
%  intermediates.tgz tar zipped directory.
%  Additionally, intermediates can be purged, greatly reducing utillized
%  space.
%  All .log .txt .mat and transformation parameter files are kept.
%  Final output volumes (.nii or .nii.gz) are kept for later visualization.
%
%  Requirements:
%       A system_and_sample_set_up.m script used to process the data, with
%       the subjectList containing IDs for which figures are to be
%       generated.
%
%  Contributors:
%           Evgeny Chumin, Indiana University, Bloomington, 2019
%%
% USER INPUT:
% system and sample set up file
sssu = '/N/dc2/projects/brainconnectomics/chumin_preproc/code/system_and_sample_set_up.m';
%
% Set to 1 if you want to purge (PERMANENTLY DELETE) the intermediates.
purge = 0;
%
%%
run(sssu);

%%
for k=1:length(subjectList)
    disp(subjectList(k).name)
    paths.subject=fullfile(paths.data,subjectList(k).name); % path to subject
    
    %% T1 clean-up
    SubjT1=fullfile(paths.subject,configs.T1.dir); % path to T1 dir
    SubjReg=fullfile(SubjT1,'registration');
    tmp=fullfile(SubjT1,'tmp');
    if ~exist(tmp,'dir')
        mkdir(tmp)
        % List of files to be kep
        source{1,1} = fullfile(SubjT1,'*log');
        source{2,1} = fullfile(SubjT1,'*txt');
        source{3,1} = fullfile(SubjT1,'*mat');
        source{4,1} = fullfile(SubjT1,'*json');
        source{5,1} = fullfile(SubjT1,'T1_fov_denoised.nii');
        source{6,1} = fullfile(SubjT1,'T1_GM_parc*.nii.gz');
        source{7,1} = fullfile(SubjT1,'T1_brain_mask_filled.nii.gz');
        source{8,1} = fullfile(SubjReg,'*mat');
        source{9,1} = fullfile(SubjReg,'*log');
        source{10,1} = fullfile(SubjReg,'T12MNI_warp.nii.gz');
        source{11,1} = fullfile(SubjReg,'MNI2T1_warp.nii.gz');
        % move kept files into tmp directory
        for s=1:length(source)
            [SUCCESS,MESSAGE,~]=movefile(source{s,1},tmp);
            if SUCCESS==0
                disp(MESSAGE)
            end
        end
        
        intermediates=fullfile(SubjT1,'intermediates');
        if ~exist(intermediates,'dir')
            mkdir(intermediates)
            % Move all remaining files into the intermediates directory    
            [SUSSESS,MESSAGE,~]=movefile(fullfile(SubjT1,'*nii*'),intermediates); %#ok<*ASGLU>
            if SUCCESS==0
                disp(MESSAGE)
            end
            [SUSSESS,MESSAGE,~]=movefile(fullfile(SubjT1,'T1_denoised.anat'),intermediates);
            if SUCCESS==0
                disp(MESSAGE)
            end
            [SUSSESS,MESSAGE,~]=movefile(fullfile(SubjT1,'registration'),intermediates);
            if SUCCESS==0
                disp(MESSAGE)
            end
            % Tar zip the intermediates directory
            archive=fullfile(SubjT1,'intermediates.tgz');
            [status,result]=system(sprintf('tar -czf %s %s',archive,intermediates));
        
            if ~isempty(result) % check for tar zip errors
                disp(result)
                disp('')
                disp('WILL NOT DELETE INTERMEDAITES')
                disp('Check that tar ran properly')
            else % remove the intermediates directory
                rmdir(fullfile(SubjT1,'intermediates'),'s')
            end
            
            % move files out of the tmp directory
            [SUCCESS,MESSAGE,~]=movefile(fullfile(tmp,'*'),SubjT1);    
            if SUCCESS==0
                disp(MESSAGE)
            end
            rmdir(tmp)
            % move transformations back into the registration directory
            mkdir(SubjReg)
            [SUCCESS,MESSAGE,~]=movefile(fullfile(SubjT1,'*dof*'),SubjReg);
            if SUCCESS==0
                disp(MESSAGE)
            end
            [SUCCESS,MESSAGE,~]=movefile(fullfile(SubjT1,'*warp*'),SubjReg);
            if SUCCESS==0
                disp(MESSAGE)
            end
        else
            disp('T1 intermediates dir exists! Check for possible errors,')
            disp('perhaps from a preceeding clean-up run.')
        end
    else
        disp('T1 tmp dir exists! Check for possible errors,')
        disp('perhaps from a preceeding clean-up run.')
    end

if purge == 1
    delete(archive)
end

    %% EPI clean-up
    
    
    
    %% DWI clean-up


end














