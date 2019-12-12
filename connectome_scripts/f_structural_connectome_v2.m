function [paths,flags,configs,parcs]=f_structural_connectome_v2(paths,flags,configs,parcs)
%                       F_STRUCTURAL_CONNECTOME
% Does the registration of anatomical and parcellation data, refrorms
% tractography, and generates connectivity matrices of sctructural data.
%
% Contributors:
%   Evgeny Chumin, Indiana University School of Medicine / IU Bloomington
% 
paths.DWI.EDDY=fullfile(paths.DWI.dir,'EDDY');
if exist(paths.DWI.EDDY,'dir')
    paths.DWI.DTIfit=fullfile(paths.DWI.dir,'DTIfit');
    if ~exist(paths.DWI.DTIfit,'dir')
        warning('Path to DTIfit directory does not exist. Exiting...')
        return
    end
else
    warning('Path to EDDY directory does not exist. Exiting...')
    return
end

%% registration of B0 to T1
if flags.DWI.regT1_2DWI==1
    disp('------------------------')
    disp('Registration of T1 to b0')
    disp('------------------------')
    
    % down sample FA image
    fileFA2mm = fullfile(paths.DWI.DTIfit,'3_DWI_FA.nii.gz');
    fileFA1mm = fullfile(paths.DWI.DTIfit,'3_DWI_FA_1mm.nii.gz');
    sentence=sprintf('%s/flirt -in %s -ref %s -out %s -applyisoxfm 1',...
        paths.FSL,fileFA2mm,fileFA2mm,fileFA1mm);
    [~,result]=system(sentence);
    
    % rigid body of T1 to b0
    disp('rigid body dof 6 to T1')
    fileIn = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    fileMat1 = fullfile(paths.DWI.dir,'T1_2_FA_dof6.mat');
    fileOut = fullfile(paths.DWI.dir,'rT1_dof6.nii.gz');
    sentence = sprintf('%s/flirt -in %s -ref %s -omat %s -dof 6 -interp spline -out %s',...
        paths.FSL,fileIn,fileFA1mm,fileMat1,fileOut);
    [~,result] = system(sentence); %#ok<*ASGLU>
    sentence=sprintf('%s/fslmaths %s -thr 0 %s',paths.FSL,fileOut,fileOut);
    [~,result]=system(sentence);
    fileMask = fullfile(paths.T1.dir,'T1_brain_mask.nii.gz');
    fileMaskDWI = fullfile(paths.DWI.dir,'rT1_brain_mask.nii.gz');
    sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
        paths.FSL,fileMask,fileFA1mm,fileMat1,fileMaskDWI);
    [~,result] = system(sentence);
    sentence=sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileOut,fileMaskDWI,fileOut);
    [~,result] = system(sentence);
    
    % aply to GM nodal parcellation images
    for k=1:length(parcs.pdir)
        if parcs.pnodal(k).true == 1    
            fileGMparc = fullfile(paths.T1.dir,sprintf('T1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
            fileOut = fullfile(paths.DWI.dir,sprintf('rT1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
            sentence=sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                        paths.FSL,fileGMparc,fileFA1mm,fileMat1,fileOut);
            [~,result] = system(sentence);
        end
    end 
end

if flags.DWI.MRtrix==1
    paths.DWI.mrtrix=fullfile(paths.DWI.dir,'MRtrix');
        if ~exist(paths.DWI.mrtrix,'dir')
            mkdir(paths.DWI.mrtrix)
        end

%% Response Function Estimation    
    disp('------------------------------')
    disp('MRtrix Streamline Tractography')
    disp('------------------------------')
    disp('1 - Estimating Response function')
    voxelsOUT=fullfile(paths.DWI.mrtrix,'csd_selected_voxels.nii.gz');
    maskIN=fullfile(paths.DWI.EDDY,'b0_brain_mask.nii.gz');
    bvecIN=fullfile(paths.DWI.EDDY,'eddy_output.eddy_rotated_bvecs');
    bvalIN=fullfile(paths.DWI.DTIfit,'3_DWI.bval');
    dataIN=fullfile(paths.DWI.EDDY,'eddy_output.nii.gz');
    % !!! name should be ties to the config below
    fileResponse=fullfile(paths.DWI.mrtrix,'tournier_response.txt');
    
    % !!! CONFIG NEEDS ADDED:
    % tournier is one of several algorith options; best choice can vary
    % depending on data quality (i.e. number of shells and directions)
    sentence1=sprintf('LD_LIBRARY_PATH= dwi2response tournier -voxels %s -force -mask %s -fslgrad %s %s %s %s',...
        voxelsOUT,maskIN,bvecIN,bvalIN,dataIN,fileResponse);
    [~,dwi2response_return]=system(sentence1);
        
%% constrained shperical deconvolution
    disp('2 - Constrained Spherical Deconvolution')
    fileFOD=fullfile(paths.DWI.mrtrix,'csd_fod.mif');
   
   % !!! CONFIG NEEDS ADDED:
   % csd is one possible algorithm for csd estimation; best one is data
   % dependent
   % additional options may be needed for multi-shell data
    sentence2=sprintf('LD_LIBRARY_PATH= dwi2fod csd -force -fslgrad %s %s -mask %s %s %s %s',...
            bvecIN,bvalIN,maskIN,dataIN,fileResponse,fileFOD);
    [~,csd_return]=system(sentence2);
        
%% ACT tissue-type volume generation
    disp('3 - Anatomically Constrained Tractography')
    brainIN=fullfile(paths.DWI.dir,'rT1_dof6.nii.gz');
    file5tt=fullfile(paths.DWI.mrtrix,'fsl5tt.nii.gz');
    
    % act needs distortion corrected data; should work with no dist corr,
    % but with nonlinear reg, but I havent tried it. 
    sentence3=sprintf('LD_LIBRARY_PATH= 5ttgen fsl -force -premasked %s %s',brainIN,file5tt);
    [~,fsl5tt_return]=system(sentence3);
        
%% generate streamlines
    disp('4 - Generating Streamlines')
    %CONFIG: 10million streamlines could be user set to other numbers
    fileStreamlines=fullfile(paths.DWI.mrtrix,'10m_streamlines.tck');
    % CONFIG: 10M can be changed
    %CONFIG: iFOD2 can be changed to other algorithms, but best one depends
    %on the data
    % CONFIG : -seed_dynamic can be switched out for other options
    % CONFIG : act can be options if data does not allow it. 
    % There may be other options withing tckgen that could be useful, this
    % is just basic usage. 
    sentence4=sprintf('LD_LIBRARY_PATH= tckgen %s %s -act %s -crop_at_gmwmi -algorithm iFOD2 -seed_dynamic %s -select 10M',...
        fileFOD,fileStreamlines,file5tt,fileFOD);
    [~,tckgen_return]=system(sentence4);
%% set min-max length boundaries
    disp('   Apply Length Filter')
    fileStreamlines2=fullfile(paths.DWI.mrtrix,'10m_10-200l_streamlines.tck');
    % CONFIG: minimum and maximum streamline lengths can be user set
    sentence5=sprintf('LD_LIBRARY_PATH= tckedit -force -minlength 10 -maxlength 200 %s %s',fileStreamlines,fileStreamlines2);
    [~,tckedit_return]=system(sentence5);
    if exist(fileStreamlines2,'file')
        delete(fileStreamlines)
    end
    
%% filter streamlines
    disp('5 - Running SIFT Filtering')
    fileFiltStreamlines=fullfile(paths.DWI.mrtrix,'1m_sift_streamlines.tck');
    % For SIFT ACT is pretty much a requirement, so if ACT cant be done, then 
    % sift shouldnt be done and tckgen can be done with less streamlines.
    sentence6=sprintf('LD_LIBRARY_PATH= tcksift -force -act %s -term_number 1M %s %s %s',...
        file5tt,fileStreamlines2,fileFOD,fileFiltStreamlines);
    [~,sift_return]=system(sentence6);    
    
end

%% Structural Connectivity Matrix Assembly
if flags.DWI.connMatrix==1
    disp('----------------------------')
    disp('Connectivity Matrix Assembly')
    disp('----------------------------')  
    fileFiltStreamlines=fullfile(paths.DWI.mrtrix,'1m_sift_streamlines.tck');
    paths.DWI.matrices=fullfile(paths.DWI.dir,'CONNmats');
        if ~exist(paths.DWI.matrices,'dir')
            mkdir(paths.DWI.matrices)
        end
    for k=1:length(parcs.pdir)
        if parcs.pnodal(k).true == 1    
            fileParc = fullfile(paths.DWI.dir,sprintf('rT1_GM_parc_%s.nii.gz',parcs.plabel(k).name));
            fileConnMatrix=fullfile(paths.DWI.matrices,sprintf('1M_2radial_density_%s.csv',parcs.plabel(k).name));
            % CONFIG assignment_radial_search can be user set
            % CONFIG scale_invnodevol: other options are available for edge
            % assignment
            % CONFIG symmetric could be optional, but its good practive to
            % have it
            % CONFIG: zero_diagonal can be optional
            sentence7{k}=sprintf('LD_LIBRARY_PATH= tck2connectome -assignment_radial_search 2 -scale_invnodevol -symmetric -zero_diagonal -force %s %s %s',...
                fileFiltStreamlines,fileParc,fileConnMatrix);
            [~,connectome_return{k}]=system(sentence7{k});
        end
    end
end
% all of the commands (sentence variables) are number postpended, so it
% will be useful to save them out into a log file. The same can be done for
% the returns from the system command. 