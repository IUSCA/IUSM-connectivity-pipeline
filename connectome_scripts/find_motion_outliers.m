function find_motion_outliers(path2uncorImage,paths)



%-------------------------------------------------------------------------%
    %% Framewise Displacement regressor
    disp('Computing fd regressors')
    
    fileOut = fullfile(paths.EPI.dir,'motionRegressor_FD.txt');
    fileMetric = fullfile(paths.EPI.dir,'motionMetric_FD.txt');
    filePlot = fullfile(paths.EPI.dir,'motionPlot_FD.png');

    if exist(fileMetric,'file')
       delete(fileMetric);
    end

    disp('Will use box-plot cutoff = P75 + 1.5*IQR')
    sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --fd',...
        paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot);

    [~,result] = system(sentence);
    disp(result)
    
%-------------------------------------------------------------------------%
    %% DVARS
    disp('Computing dvars regressors')
    
    fileOut = fullfile(paths.EPI.dir,'motionRegressor_DVARS.txt');
    fileMetric = fullfile(paths.EPI.dir,'motionMetric_DVARS.txt');
    filePlot = fullfile(paths.EPI.dir,'motionPlot_DVARS.png');
    if exist(fileMetric,'file')
        delete(fileMetric);
    end
    
    disp('Will use box-plot cutoff = P75 + 1.5*IQR')
    sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --dvars',...
        paths.FSL,path2uncorImage,fileOut,fileMetric,filePlot);

    [~,result] = system(sentence);
    disp(result)
end
