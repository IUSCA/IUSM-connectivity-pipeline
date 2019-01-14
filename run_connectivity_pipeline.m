function run_connectivity_pipeline
%                       RUN_CONNECTIVITY_PIPELINE
% Executes processing of anatomical, functional, and diffusion MRI date as
% defined in the set_up scripts:
%       SYSTEM_AND_SAMPLE_SET_UP.m
%       BATCH_SET_UP.m 
%
% The two scripts are identified via matlab gui and ran, thus initializing
% the environment for the pipeline.
%
%      Evgeny Chumin, Indiana University School of Medicine, 2018
%      John West, Indiana University School of Medicine, 2018
%     Mario Dzemidzic, Indiana University School of Medicine, 2018
%       Zikai Lin, Indiana University School of Medicine, 2018
%
%%
% 
% JDW edit 03/16/2018 - adjusted get_readout to use user defined DICOM extension
%
%%
% Set path/name to the batch set_up script
[filename_b,pathname_b]=uigetfile('*.m','Select Batch Setup .m File');
run(fullfile(pathname_b,filename_b));
% Set path/name to the system and sample set_up script
[filename_s,pathname_s]=uigetfile('*.m','Select System/Sample Setup .m File');
run(fullfile(pathname_s,filename_s));
% save name in paths structure for reference
paths.batch=fullfile(pathname_b,filename_b);
% save name in paths structure for reference
paths.setup=fullfile(pathname_s,filename_s);

% Export FSL
[~,~]=system(FSLsetup); clear FSLsetup


%%
if flags.global.T1_prepare_A==1
    % Run T1_A on all subjects in subjectList
for i=1:length(subjectList)
    run(fullfile(pathname_b,filename_b)) % initialize default configs
    paths.subject = fullfile(paths.data,subjectList(i).name);
    % check that subject path and T1 directory exist
    if exist(paths.subject,'dir') && exist(fullfile(paths.subject,configs.name.T1),'dir') %#ok<*NODEF>
        paths.T1.dir = fullfile(paths.subject,configs.name.T1);
        diary(fullfile(paths.T1.dir,sprintf('diary_%s.log',datestr(now,'yyyymmdd'))))
        % check that T1 contains something (whether its dicoms dir or a nifti)
        if ~isempty(paths.T1.dir)
            disp('------------------------')
            fprintf('T1_prepare_A on %s\n',subjectList(i).name)
            disp('------------------------')
            % run preprocessing
            [paths,flags,configs]=f_T1_prepare_A(paths,flags,configs);
            % save the cofiguration variables for this run
            configFile=fullfile(paths.subject,sprintf('configs_T1_A_%s.mat',datestr(now,'yyyymmdd')));
            if exist(configFile,'file')
                count=dir(fullfile(paths.subject,sprintf('configs_T1_A_%s*',datestr(now,'yyyymmdd'))));
                count=length(count);
                configFile=fullfile(paths.subject,sprintf('configs_T1_A_%s_run%d.mat',datestr(now,'yyyymmdd'),count+1));
            end
            save(configFile,'-struct','configs','T1');
        else
            warning('T1 directory is empty; skipping subject')
            return
        end
    else
        warning('Either subject or contained T1 directory does not exit; skipping subject')
        return
    end
    diary off
end
end

%%
if flags.global.T1_prepare_B==1
    % Run T1_B on all subjects in subjectList
for i=1:length(subjectList)
    run(fullfile(pathname_b,filename_b)) % initialize default configs
    paths.subject = fullfile(paths.data,subjectList(i).name);
    % check that subject path and T1 directory exist
    paths.T1.dir = fullfile(paths.subject,configs.name.T1);
    if exist(paths.subject,'dir') && exist(paths.T1.dir,'dir')
        diary(fullfile(paths.T1.dir,sprintf('diary_%s.log',datestr(now,'yyyymmdd'))))
        disp('------------------------')
        fprintf('T1_prepare_B on %s\n',subjectList(i).name)
        disp('------------------------')
        % run registration and segmentation
        [paths,flags,configs,parcs]=f_T1_prepare_B(paths,flags,configs,parcs);
        % save the cofiguration variables for this run
        configFile=fullfile(paths.subject,sprintf('configs_T1_B_%s.mat',datestr(now,'yyyymmdd')));
        if exist(configFile,'file')
            count=dir(fullfile(paths.subject,sprintf('configs_T1_B_%s*',datestr(now,'yyyymmdd'))));
            count=length(count);
            configFile=fullfile(paths.subject,sprintf('configs_T1_B_%s_run%d.mat',datestr(now,'yyyymmdd'),count+1));
        end
        save(configFile,'-struct','configs','T1');
    else
        warning('Either subject or contained T1 directory does not exit; skipping subject')
        return
    end
    diary off
end
end

%%
if flags.global.fMRI_A==1
   % run EPI processing on all subjects in subjectList
for i=1:length(subjectList)
    run(fullfile(pathname_b,filename_b)) % initialize default configs
    paths.subject = fullfile(paths.data,subjectList(i).name);
    paths.T1.dir=fullfile(paths.subject,configs.name.T1);
    % check that subject path and T1 directory exist
    if exist(paths.subject,'dir') && exist(paths.T1.dir,'dir')
        % Generate list of EPI scan directories
        dircont=dir(paths.subject); dircont(1:2)=[];% find content of subject dir
        epiList=struct.empty;
        for e=1:length(dircont)
            if dircont(e).isdir==1 && ~isempty(strfind(dircont(e).name,configs.name.epiFolder))
                epiList(end+1).name=dircont(e).name; %#ok<*AGROW>
            end
        end
        if ~isempty(epiList)
            configs.EPI.sessions=epiList.name;
            % For each session under this subject    
            for j=1:length(epiList)
                run(fullfile(pathname_b,filename_b)) % initialize default configs
            % Operating on the scans set in configs
                if j >= configs.EPI.epiMin && j <= configs.EPI.epiMax
                    paths.EPI.dir=fullfile(paths.subject,epiList(j).name);
                    diary(fullfile(paths.EPI.dir,strcat('diary_',datestr(now,'yyyymmdd'),'.log')))
                    disp('---------------------------------')
                    fprintf('fMRI_A on %s\n',subjectList(i).name)
                    fprintf('EPI series: %s\n',epiList(j).name)
                    disp('---------------------------------')
                    % run fMRI processing
                    [paths,flags,configs,parcs]=f_functional_connectivity(paths,flags,configs,parcs);
                    % save the cofiguration variables for this run
                    configFile=fullfile(paths.subject,sprintf('configs_%s_A_%s.mat',epiList(j).name,datestr(now,'yyyymmdd')));
                    if exist(configFile,'file')
                        count=dir(fullfile(paths.subject,sprintf('configs_%s_A_%s*',epiList(j).name,datestr(now,'yyyymmdd'))));
                        count=length(count);
                        configFile=fullfile(paths.subject,sprintf('configs_%s_A_%s_run%d.mat',epiList(j).name,datestr(now,'yyyymmdd'),count+1));
                    end
                    save(configFile,'-struct','configs','EPI');
                end
            end
        else
            warning('epiList is empty or does not exit. Check consistency of naming convention.')
            return
        end
    else
        warning('Either subject or contained T1 directory does not exit; skipping subject')
        return
    end
    diary off
end 
end

%%
if flags.global.fMRI_B==1   
    % Generate figures for all subjects in subjectList
    for i=1:length(subjectList)
        run(fullfile(pathname_b,filename_b)) % initialize default configs
        paths.subject = fullfile(paths.data,subjectList(i).name);
        paths.T1.dir=fullfile(paths.subject,configs.name.T1);
        % check that subject path and T1 directory exist
        if exist(paths.subject,'dir') && exist(paths.T1.dir,'dir')
            % Generate list of EPI scan directories
            dircont=dir(paths.subject); dircont(1:2)=[];% find content of subject dir
            epiList=struct.empty;
            for e=1:length(dircont)
                if dircont(e).isdir==1 && ~isempty(strfind(dircont(e).name,configs.name.epiFolder))
                    epiList(end+1).name=dircont(e).name; %#ok<*AGROW>
                end
            end
            if ~isempty(epiList)
                configs.EPI.sessions=epiList;
                % For each session under this subject    
                for j=1:length(epiList)
                    % Operating on the scans set in configs
                    if j >= configs.EPI.epiMin && j <= configs.EPI.epiMax
                        paths.EPI.dir=fullfile(paths.subject,epiList(j).name);
                        diary(fullfile(paths.EPI.dir,strcat('diary_',datestr(now,'yyyymmdd'),'.log')))
                        if flags.EPI.GS==1
                            paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
                            GSreg_name = 'GSreg\_yes';
                        elseif flags.EPI.GS==0
                            paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
                            GSreg_name = 'GSreg\_no';
                        else
                            warning('flags.EPI.GS not specified. Exiting...')
                        end
                        if exist(paths.EPI.epiGS,'dir')
                            subjectinfo = strcat(subjectList(i).name,', ',epiList(j).name,', ',GSreg_name);
                            disp('------------------------')
                            fprintf('fMRI_B on %s\n',subjectinfo)
                            fprintf('EPI series: %s\n',epiList(j).name)
                            disp('------------------------')
                            % Generate the figures
                            [paths,flags,configs,parcs]=f_evaluateFC(paths,flags,configs,parcs,subjectinfo);
                            close all;
                            % save the cofiguration variables for this run
                            configFile=fullfile(paths.subject,sprintf('configs_%s_B_%s.mat',epiList(j).name,datestr(now,'yyyymmdd')));
                            if exist(configFile,'file')
                                count=dir(fullfile(paths.subject,sprintf('configs_%s_B_%s*',epiList(j).name,datestr(now,'yyyymmdd'))));
                                count=length(count);
                                configFile=fullfile(paths.subject,sprintf('configs_%s_B_%s_run%d.mat',epiList(j).name,datestr(now,'yyyymmdd'),count+1));
                            end
                            save(configFile,'-struct','configs','EPI');
                        else
                            warning('Either GSReg_yes or _no directory does not exit; skipping subject')
                            return
                        end
                        diary off
                    end
                end
            else
                warning('No %s containing session names found. Check for consistency of nomenclature',configs.name.epiFolder)
            end
        else
            warning('Either path to subject or %s not found',configs.name.epiFolder)
        end
    end
end

%%
if flags.global.DWI_A==1
    for i=1:length(subjectList) % For each subject
        run(fullfile(pathname_b,filename_b)) % initialize default configs
        % Set up environment
        paths.subject = fullfile(paths.data,subjectList(i).name);
        paths.T1.dir = fullfile(paths.subject,configs.name.T1);
        paths.DWI.dir = fullfile(paths.subject,configs.name.DWI);
        if exist(paths.DWI.dir,'dir')
            diary(fullfile(paths.DWI.dir,strcat('diary_',datestr(now,'yyyymmdd'),'.log')))
            disp('---------------------------')
            fprintf('Processing DWI of %s\n',subjectList(i).name)
            disp('---------------------------')
            if isempty(configs.DWI.readout)
            % calculate readout time
            [configs.DWI.readout]=get_readout(paths,configs.name.dcmFiles);
            end
            % run DWI preprocessing
            [paths,flags,configs]=f_preproc_DWI(paths,flags,configs);
            % save the cofiguration variables for this run
            configFile=fullfile(paths.subject,sprintf('configs_DWI_A_%s.mat',datestr(now,'yyyymmdd')));
                if exist(configFile,'file')
                    count=dir(fullfile(paths.subject,sprintf('configs_DWI_A_%s*',datestr(now,'yyyymmdd'))));
                    count=length(count);
                    configFile=fullfile(paths.subject,sprintf('configs_DWI_A_%s_run%d.mat',datestr(now,'yyyymmdd'),count+1));
                end
            save(configFile,'-struct','configs','DWI');
        else
            warning('Subject DWI directory does not exit; skipping subject')
            return
        end
        diary off
    end
end

%%
if flags.global.DWI_B==1
    for i=1:length(subjectList) % For each subject
        run(fullfile(pathname_b,filename_b)) % initialize default configs
        % Set up environment
        paths.subject = fullfile(paths.data,subjectList(i).name);
        paths.T1.dir = fullfile(paths.subject,configs.name.T1);
        paths.DWI.dir = fullfile(paths.subject,configs.name.DWI);
        if exist(paths.DWI.dir,'dir')
            diary(fullfile(paths.DWI.dir,strcat('diary_',datestr(now,'yyyymmdd'),'.log')))
            disp('---------------------------')
            fprintf('Processing DWI of %s\n',subjectList(i).name)
            disp('---------------------------')
            % run DWI preprocessing
            [paths,flags,configs,parcs]=f_structural_connectome(paths,flags,configs,parcs);
            % save the cofiguration variables for this run
            configFile=fullfile(paths.subject,sprintf('configs_DWI_B_%s.mat',datestr(now,'yyyymmdd')));
            if exist(configFile,'file')
                    count=dir(fullfile(paths.subject,sprintf('configs_DWI_B_%s*',datestr(now,'yyyymmdd'))));
                    count=length(count);
                    configFile=fullfile(paths.subject,sprintf('configs_DWI_B_%s_run%d.mat',datestr(now,'yyyymmdd'),count+1));
            end
            save(configFile,'-struct','configs','DWI');
        else
            warning('Subject DWI directory does not exit; skipping subject')
            return
        end
        diary off
    end
end

%%
if flags.global.DWI_C==1
    for i=1:length(subjectList) % For each subject
        run(fullfile(pathname_b,filename_b)) % initialize default configs
        % Set up environment
        paths.subject = fullfile(paths.data,subjectList(i).name);
        paths.T1.dir = fullfile(paths.subject,configs.name.T1);
        paths.DWI.dir = fullfile(paths.subject,configs.name.DWI);
        if exist(paths.DWI.dir,'dir')
            diary(fullfile(paths.DWI.dir,strcat('diary_',datestr(now,'yyyymmdd'),'.log')))
            disp('---------------------------')
            fprintf('Processing DWI_C of %s\n',subjectList(i).name)
            disp('---------------------------')
            % run probabilistic tractography
            [paths,flags,configs,parcs]=f_prob_conn(paths,flags,configs,parcs,params);
            % save the cofiguration variables for this run
            configFile=fullfile(paths.subject,sprintf('configs_DWI_C_%s.mat',datestr(now,'yyyymmdd')));
            if exist(configFile,'file')
                    count=dir(fullfile(paths.subject,sprintf('configs_DWI_C_%s*',datestr(now,'yyyymmdd'))));
                    count=length(count);
                    configFile=fullfile(paths.subject,sprintf('configs_DWI_C_%s_run%d.mat',datestr(now,'yyyymmdd'),count+1));
            end
            save(configFile,'-struct','configs','DWI');
        else
            warning('Subject DWI directory does not exit; skipping subject')
            return
        end
        diary off
    end
end