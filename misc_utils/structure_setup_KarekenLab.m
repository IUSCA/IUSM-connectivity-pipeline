%% Structural Data Configuration

paths.datadir = '/bgodata01/BAD2/CONNECTIVITY/datadir_091418/'; % target directory, input directory for pipeline
paths.basePath = '/bgodata01/BAD2/';
paths.source_path = '/bgodata01/BAD2/'; % Path to subject source data

configs.removeExistDir = 1;
% "BAD0159","BAD0160", "BAD0171" excluded for now, already set up
subjectList = ["BAD0182",...
    "BAD0189", "BAD0196", "BAD0197", "BAD0201", "BAD0203", "BAD0204",...
    "BAD0215", "BAD0217", "BAD0219", "BAD0220", "BAD0221",...
    "BAD0231", "BAD0235", "BAD0249", "BAD0256", "BAD0257", "BAD0270", "BAD0277",...
    "BAD0280", "BAD0286", "BAD0312", "BAD0320", "BAD0326"]; % A list of Subject ID 


params.nsub = length(subjectList); % Number of subjects
params.nEPI = 5; % Total number of task includes coonn



taskList = ["CONN1", "TASK1", "TASK2", "TASK3", "TASK4"]; % Depends on the design of the study
diary('structure_set_up_diary_20180919'); 
% Remove the datadir that existed in target directory
if configs.removeExistDir == 1
    for sub = 1:length(subjectList)
        sentence = sprintf("rm -rd %s%s", paths.datadir, subjectList(sub));
        system(sentence);
    end
end

%% Setup datadir Structure
if exist(paths.datadir, 'dir') == 7
    disp('datadir exist, please manully remove it, with cautious.')
else 
    disp('datadir not exists, setting up structure...')
    mkdir(paths.datadir);
end

% Loop over subject
for subi = 1:params.nsub
    paths.SubjectSrcPath = fullfile(paths.basePath, subjectList(subi));
    sentence = sprintf("cd %s", paths.SubjectSrcPath);
    system(sentence);
    
    % Check how many tasks does this subject conducted.
    params.nEPI = 0;
    params.subEPItask = [];
    for task = 1:length(taskList)
        if exist(strcat(paths.SubjectSrcPath,"/MBEPI-",taskList(task),"/"), 'dir')
            params.nEPI = params.nEPI + 1; % Found 1 more task
            params.subEPItask = [params.subEPItask, taskList(task)];
        else
            fprintf("For subject %s, %s is missing\n", subjectList(subi), taskList(task));
        end
    end
    
    if ~ismember("CONN1", params.subEPItask)
        error(sprintf("CONN1 is not in subject %s, please fix it all exclude it from the subjectList", subjectList(subi)));
    end
    
    
    % T1 structure set up
    paths.T1.dir = fullfile(paths.datadir, subjectList(subi),"T1/");
    paths.T1.dcm = fullfile(paths.T1.dir, "DICOMS");
    mkdir(paths.T1.dir);
    mkdir(paths.T1.dcm);
    if exist(paths.T1.dir, 'dir') == 7
        disp('T1 directory created.')
    else
        warning('Cannot make T1 directory.')
    end
    
    
    
    paths.EPI.epi_paths = {}; % should contains at least nepi directory
    paths.EPI.epi_dcm = {};
    paths.EPI.epi_unwarp = {};
    paths.EPI.AP = {};
    paths.EPI.PA = {};
 

    % EPI structure set up
    for i = 1:params.nEPI
        epi_n = strcat("EPI", num2str(i));
        epi_dir = fullfile(paths.datadir, subjectList(subi), epi_n);
        paths.EPI.epi_paths{i} = epi_dir;
        paths.EPI.epi_dcm{i} = fullfile(paths.EPI.epi_paths{i}, "DICOMS");
        paths.EPI.epi_unwarp{i} = fullfile(paths.EPI.epi_paths{i}, "UNWARP");
        paths.EPI.AP{i} = fullfile(paths.EPI.epi_unwarp{i}, "SEFM_AP_DICOMS");
        paths.EPI.PA{i} = fullfile(paths.EPI.epi_unwarp{i}, "SEFM_PA_DICOMS");
        
        
        mkdir(epi_dir);
        mkdir(paths.EPI.epi_dcm{i});
        mkdir(paths.EPI.epi_unwarp{i});
        mkdir(paths.EPI.AP{i});
        mkdir(paths.EPI.PA{i});
        
        if exist(paths.EPI.epi_paths{i}, 'dir') == 7
            mssg = sprintf("EPI%d directory created.", i);
            disp(mssg)
        else
            mssg = sprintf("Error while creating EPI%d directory.", i);
            warning(mssg)
        end 
    end
        
    % DWI structure set up
    paths.DWI.dir = fullfile(paths.datadir,subjectList(subi), "DWI");
    paths.DWI.dcm = fullfile(paths.DWI.dir, "DICOMS");
    paths.DWI.unwarp = fullfile(paths.DWI.dir, "UNWARP");
    paths.DWI.B0_PA = fullfile(paths.DWI.unwarp, "B0_PA_DCM");
    mkdir(paths.DWI.dir);
    mkdir(paths.DWI.dcm);
    mkdir(paths.DWI.unwarp);
    mkdir(paths.DWI.B0_PA);
    if exist(paths.DWI.dir, 'dir') == 7
        disp('DWI directory created.')
    else
        warning('Cannot make T1 directory.')
    end      


    
    %% copy dicoms files to T1
    disp("==================================")
    disp("Setting up T1 ... ")
    sentence = sprintf("cp %s/MPRAGE1/*.dcm %s", paths.SubjectSrcPath, paths.T1.dcm);
    [status, mssg] = system(sentence);
    if status == 0
        disp("T1 DICOMS copied")
    else
        warning(mssg);
    end

    
    %% copy dicoms files to EPI
    for i = 1:params.nEPI
        epi_dir = paths.EPI.epi_paths{i};
        
        
        if i == 1
            % First set up CONN1
            disp("==================================")
            disp("setting up EPI CONN1...")
            sentence = sprintf("cp %s/MBEPI-CONN1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.epi_dcm{i});
            [status, mssg] = system(sentence);
            if status == 0
                disp_mssg = sprintf("EPI%d DICOMS files copied.", i);
                disp(disp_mssg)
            else
                disp_mssg = sprintf(mssg);
                warning(disp_mssg);
            end 

            % Copy FIELD-MAP AP
            sentence = sprintf("cp %s/SEFIELD-AP-CONN1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.AP{i});
            [status, mssg] = system(sentence);
            if status == 0
                disp("EPI CONN1 SPEFIELD AP MAP copied")
            else
                warning(mssg);
            end
            
            % Copy FIELD-MAP PA
            sentence = sprintf("cp %s/SEFIELD-PA-CONN1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.PA{i});
            [status, mssg] = system(sentence);
            if status == 0
                disp("EPI CONN1 SPEFIELD PA MAP copied")
            else
                warning(mssg);
            end 
        else
            
           % Setting up EPI-task
           disp("==================================")
           fprintf("setting up EPI task %d...\n", i-1)
           
           % Task DICOMS
           
            sentence = sprintf("cp %s/MBEPI-%s/*.dcm %s",paths.SubjectSrcPath,params.subEPItask(i), paths.EPI.epi_dcm{i});
            [status, mssg] = system(sentence);
            if status == 0
                disp_mssg = sprintf("EPI%d DICOMS files copied.", i);
                disp(disp_mssg)
            else
                warning(mssg);
            end 

            % Copy FIELD-MAP
            
            % CONN1 or TASK1?
            sentence = sprintf("readlink -f %s/MBEPI-%s",paths.SubjectSrcPath, params.subEPItask(i)); % readlink from task id
            [status, link] = system(sentence);
            link = strsplit(link, '/');
            series = link{5}(2:(length(link{5})-1)); % The last char is \n
            
            if (4 < str2num(series)) &&  (str2num(series) < 16)
                % FIELD-MAP AP CONN
                sentence = sprintf("cp %s/SEFIELD-AP-CONN1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.AP{i});
                [status, mssg] = system(sentence);
                if status == 0
                    fprintf("Series %s EPI CONN1 SPEFIELD AP MAP copied to EPI%d\n", series, i);
                else
                    warning(mssg);
                end
                
               % FIELD-MAP PA CONN
                sentence = sprintf("cp %s/SEFIELD-PA-CONN1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.PA{i});
                [status, mssg] = system(sentence);
                if status == 0
                    fprintf("Series %s EPI CONN1 SPEFIELD PA MAP copied to EPI%d\n", series, i);
                else
                    warning(mssg);
                end   
                
            elseif str2num(series) > 16
                % FIELD-MAP AP TASK
                sentence = sprintf("cp %s/SEFIELD-AP-TASK1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.AP{i});
                [status, mssg] = system(sentence);
                if status == 0
                    fprintf("Series %s EPI TASK1 SPEFIELD AP MAP copied to EPI%d\n", series, i);
                else
                    warning(mssg);    
                end 
                
                % FIELD-MAP PA TASK
                sentence = sprintf("cp %s/SEFIELD-PA-TASK1/*.dcm %s",paths.SubjectSrcPath, paths.EPI.PA{i});
                [status, mssg] = system(sentence);
                if status == 0
                    fprintf("Series %s EPI TASK1 SPEFIELD PA MAP copied to EPI%d\n", series, i);
                else
                    warning(mssg);
                end     
            end
        end
    end
        
    disp("EPI copied!!")
    disp("================================")
    disp("Now copying DWI...");
    
    %% copy dicoms files to DWI    
    sentence = sprintf("cp %s/DTI1/*.dcm %s",paths.SubjectSrcPath, paths.DWI.dcm);
    [status, mssg] = system(sentence);
    if status == 0
        disp("DWI DICOMS copied")
    else
        warning(mssg);
    end
    
    %% copy unwarp dicom files to DWI    
    sentence = sprintf("cp %s/DTIb01/*.dcm %s",paths.SubjectSrcPath, paths.DWI.B0_PA);
    [status, mssg] = system(sentence);
    if status == 0
        disp("DWI unwarp dcm copied.")
    else
        warning(mssg);
    end  
    disp("================================")
fprintf("Finish copying data for subject %s\n", subjectList(subi))
disp("================================")
end
diary OFF;

    
    
    
    
    









