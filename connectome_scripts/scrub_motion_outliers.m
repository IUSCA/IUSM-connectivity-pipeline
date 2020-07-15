function [good_vols] = scrub_motion_outliers(paths,configs)

for t = 1:2 % 1-FD; 2-DVARS
switch t
    case 1
        if isempty(configs.EPI.FDcut)
            fileIn = fullfile(paths.EPI.dir,'motionRegressor_FD.txt');
            if exist(fileIn,'file') 
                fd_scrub = sum(load(fileIn),2); % motionMetric_fd
                fprintf('-- identified %d FD outliers based on FSL outlier criteria\n',nnz(fd_scrub))
            else
                error('file %s not found!',fileIn)
            end
        else
            fileIn = fullfile(paths.EPI.dir,'motionMetric_FD.txt');
            fd_scrub = load(fileIn) > configs.EPI.FDcut;
            fprintf('-- identified %d FD outliers based exceeding FD of %d\n',nnz(fd_scrub),configs.EPI.FDcut)
        end
    case 2
        if isempty(configs.EPI.DVARScut)
            fileIn = fullfile(paths.EPI.dir,'motionRegressor_DVARS.txt');
            if exist(fileIn,'file')
                dvars_scrub = sum(load(fileIn),2);
                fprintf('-- identified %d DVARS outliers based on FSL outlier criteria\n',nnz(dvars_scrub))       
            else
                error('file %s not found!',fileIn)
            end
        else
            fileIn = fullfile(paths.EPI.dir,'motionMetric_DVARS.txt');
            dvars_scrub = load(fileIn) > configs.EPI.DVARScut;
            fprintf('-- identified %d DVARS outliers based exceeding DVARS of %d\n',nnz(fd_scrub),configs.EPI.DVARScut)
        end
end        
end  
good_vols=~sum(horzcat(fd_scrub,dvars_scrub),2);
fprintf('-- %d total outliers to be scrubbed from dataset\n',sum(good_vols==0))
save(fullfile(paths.EPI.dir,'scrubbing_goodvols.mat'),'good_vols')
end
    

