function [fd_scrub,dvars_scrub] = find_motion_outliers(path2uncorImage,path2mask,paths,configs)



%-------------------------------------------------------------------------%
    %% Framewise Displacement regressor
    disp('Computing fd regressors')
    fileOut = fullfile(paths.EPI.dir,'motionRegressor_FD.txt');
    fileMetric = fullfile(paths.EPI.dir,'motionMetric_FD.txt');
    filePlot = fullfile(paths.EPI.dir,'motionPlot_FD.png');
    if exist(fileOut,'file')
        delete(fileOut);
    end
    if exist(fileMetric,'file')
       delete(fileMetric);
    end
    if isempty(configs.EPI.FDcut,'var')
        disp('Will use box-plot cutoff = P75 + 1.5*IQR')
        sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --fd -m %s',...
            paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot,path2mask);
    else
        sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --fd --thresh=%0.2f -m %s',...
            paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot,configs.EPI.FDcut,path2mask);
    end
    [~,result] = system(sentence);
    disp(result)
    
    if exist(fileOut,'file') 
        fd_scrub = sum(load(fileOut),2); % motionMetric_fd
        fprintf('-- identified %d FD outliers\n',nnz(fd_scrub))
    else
        error('file %s not found!',fileOut)
    end
%-------------------------------------------------------------------------%
    %% DVARS
    disp('Computing dvars regressors')

    fileOut = fullfile(paths.EPI.dir,'motionRegressor_DVARS.txt');
    fileMetric = fullfile(paths.EPI.dir,'motionMetric_DVARS.txt');
    filePlot = fullfile(paths.EPI.dir,'motionPlot_DVARS.png');
    if exist(fileMetric,'file')
        delete(fileMetric);
    end
    
    if isempty(configs.EPI.DVARScut,'var')
        disp('Will use box-plot cutoff = P75 + 1.5*IQR')
        sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --dvars -m %s',...
            paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot,path2mask);
    else
        sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --dvars --thresh=%0.2f -m %s',...
            paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot,configs.EPI.DVARScut,path2mask);
    end
    [~,result] = system(sentence);
    disp(result)
    
    if exist(fileMetric,'file')
        dvars_scrub = sum(load(fileOut),2);
        fprintf('-- identified %d DVARS outliers\n',nnz(dvars_scrub))       
    else
        error('file %s not found!',fileMetric)
    end
