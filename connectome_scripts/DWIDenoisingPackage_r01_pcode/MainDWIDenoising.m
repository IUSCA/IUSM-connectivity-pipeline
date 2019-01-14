% Pierrick Coupe - pierrick.coupe@labri.fr
% Jose V. Manjon - jmanjon@fis.upv.es
% LaBRI UMR 5800
% Universite de Bordeaux 1
%
% Copyright (C) 2008-2013 Pierrick Coupe and Jose V. Manjon
% ***************************************************************************/

warning off;
clc;
clf;
clear all;
close all;


% Number of threads,  automatic by default but can be defined by the user.
% We suggest to x2 'maxNumCompThreads' on intel architecture to get the 
% number of logical CPUs
nbthread = maxNumCompThreads*2; % INTEL with old matlab
%nbthread=feature('numcores')*2; % INTEL with matlab 2013b
%nbthread=feature('numcores'); % AMD


%Add path
addpath 'spm8';
addpath 'gui';
packagepath = pwd;
addpath(genpath(fullfile(packagepath, 'DWIDenoisingPackage')));



fprintf('%s \n', 'Welcome to DWI Denoising software:')
fprintf('%s \n\n', 'Copyright (C) 2008-2013 Pierrick Coupe and Jose V. Manjon')

[namein, pathin, filterindex] = uigetfile({  '*.nii','NIFTI image (*.nii)'}, 'Select one or several files (using +CTRL or +SHIFT)','MultiSelect', 'on');
if isequal(namein,0) | isequal(pathin,0)
    disp('User pressed cancel')
else
    
    
    [flag method beta searchradius suffixfile rician verbose] = gui;
    

    
    fprintf('%s\n\n', 'Selected parameters')
    if (method==1)
        com = sprintf('Optimized NLM Filter\n');
        disp(com)
    end
    if (method==2)
        com = sprintf('Adaptative ONLM Filter\n');
        disp(com)
    end
    if (method==3)
        com = sprintf('Multi-resolution ONLM Filter\n');
        disp(com)
    end
    
    if (method==4)
        com = sprintf('Oracle DCT Filter\n');
        disp(com)
    end
    
    if (method==5)
        com = sprintf('Prefiltered rationally invariant NLM Filter\n');
        disp(com)
    end

  if (method==6)
        com = sprintf('Local PCA filter\n');
        disp(com)
  end
    
 if((method~=6))
    patchradius = 1;
    com = sprintf('beta: %0.1f', beta);
    disp(com)
    com = sprintf('patch size: %dx%dx%d voxels', 2*patchradius+1, 2*patchradius+1, 2*patchradius+1);
    disp(com)
end
    if (rician==1)
        com = sprintf('Rician noise model\n');
        disp(com)
    else
        com = sprintf('Gaussian noise model\n');
        disp(com)
    end
    
      if (verbose==1)
        com = sprintf('Display option on\n');
        disp(com)
    else
        com = sprintf('Display option off\n');
        disp(com)
    end
    
    
    if (flag==1)
        
        if (iscell(namein))
            nbfile = size(namein,2);
        else
            nbfile = 1;
        end
        
        for f = 1 : nbfile
            
            
            
            if(nbfile>1)
                filenamein = namein{f};
            else
                filenamein = namein;
            end
            
            disp(['Input file : ', fullfile(pathin, filenamein)])
            [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
            nout=[name_s suffixfile ext];
            pathout = pathin;
            
            disp(['Output file: ', fullfile(pathout, nout)])
            
            
            
            cname = fullfile(pathin, filenamein);

            VI=spm_vol(cname);
            ima=spm_read_vols(VI);
            s=size(ima);
            ima = double(ima);

		dirname = sprintf('%s.bvec',cname(1:end-4));
            
            if(exist(dirname))
                dir = load(dirname);
                if (size(dir,2)~=3)
                    dir = dir';
                end
                if (size(dir,2)~=3)
                    disp('your format for gradient direction is not correct')
                    return
                end
            else
                
                dir = directiondetection(ima);
                if (size(dir,2)~=3)
                    dir = dir';
                end
                if (size(dir,2)~=3)
                    disp('your format for gradient direction is not correct')
                    return
                end
            end
            


            
            % fixed range
            
            map = isnan(ima(:));
            ima(map) = 0;
            map = isinf(ima(:));
            ima(map) = 0;
            mini = min(ima(:));
            ima = ima - mini;
            maxi=max(ima(:));
            ima=ima*255/maxi;
            
            
            %Noise estimation
            if((method~=2) & (method~=6))     
                [hfinal, hvect, nbB0, hobj, hbg] = DWINoiseEstimation(ima, dir, rician, verbose);
                if(isnan(hfinal))
                    disp('error during noise estimation')
                    continue;
                end
            end
            
            
            
            
            % Denoising
  
            if(method==1)
                DWIdenoised = DWIDenoisingORNLM(ima, hfinal, beta, patchradius, 5, rician, nbthread, verbose);
            end
            
            if(method==2)
                
                DWIdenoised = DWIDenoisingAONLM(ima,patchradius, 3, beta , rician, nbthread, verbose);
            end
            
            if(method==3)
                
                DWIdenoised = DWIDenoisingORNLMMultires(ima, dir, hfinal, beta, patchradius, 3, rician, nbthread, verbose);
            end
            if(method==4)
                DWIdenoised = DWIDenoisingODCT(ima, hfinal, beta, rician, nbthread, verbose);
            end
            
            if(method==5)
                DWIdenoised = DWIDenoisingPRINLM(ima, hfinal, beta, rician, nbthread, verbose);
            end
            
            if(method==6)
               
                beta
                rician
                nbthread
                verbose
                DWIdenoised = DWIDenoisingLPCA(ima, beta, rician, nbthread, verbose);
            end            
  

            
            % Original intensity range
            DWIdenoised=DWIdenoised*maxi/255 + mini;
            map = find(DWIdenoised<0);
            DWIdenoised(map) =0;
            map = isnan(DWIdenoised);
            DWIdenoised(map)=0;
            
            
            VO = VI; % copy input info for output image
            
            outfilename =fullfile(pathout, nout);
           
  	   for i=1:s(4)
                VO(i).dim=s(1:3);
                VO(i).fname = outfilename;
                spm_write_vol(VO(i),DWIdenoised(:,:,:,i));
       end
            
        disp('Filtering done')
            
        end
    end
end


