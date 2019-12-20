

load('subjecLists_skyra_prisma_con_aud.mat')

for i=1:4
    if i==1
        sample=SKYRA_CON;
    elseif i==2
        sample=SKYRA_AUD;
    elseif i==3
        sample=PRISMA_CON;
    elseif i==4
        sample=PRISMA_AUD;
    end
    
    for j=1:length(sample)
        fprintf('Running %s\n',sample(j).name)
        DWIdir=fullfile(sample(j).folder,sample(j).name,'DWI');
        T1dir=fullfile(sample(j).folder,sample(j).name,'T1');
        OUTdir=fullfile(DWIdir,'MRtrix');
        if ~exist(OUTdir,'dir')
            mkdir(OUTdir)
        end
    % response function estimation
  %  disp('1 - Estimating Response function')
        voxelsOUT=fullfile(OUTdir,'csd_selected_voxels.nii.gz');
        maskIN=fullfile(DWIdir,'EDDY/b0_brain_mask.nii.gz');
        bvecIN=fullfile(DWIdir,'EDDY/eddy_output.eddy_rotated_bvecs');
        bvalIN=fullfile(DWIdir,'DTIfit/3_DWI.bval');
        dataIN=fullfile(DWIdir,'EDDY/eddy_output.nii.gz');
        fileResponse=fullfile(OUTdir,'tournier_response.txt');
        
        sentence1=sprintf('dwi2response tournier -voxels %s -force -mask %s -fslgrad %s %s %s %s',...
            voxelsOUT,maskIN,bvecIN,bvalIN,dataIN,fileResponse);
        %[~,dwi2response_return]=system(sentence);
        
    % constrained shperical deconvolution
 %   disp('2 - Constrained Spherical Deconvolution')
        fileFOD=fullfile(OUTdir,'csd_fod.mif');
        
        sentence2=sprintf('dwi2fod csd -fslgrad %s %s -mask %s %s %s %s',...
            bvecIN,bvalIN,maskIN,dataIN,fileResponse,fileFOD);
        %[~,csd_return]=system(sentence);
        
    % ACT tissue-type volume generation
  %  disp('3 - Anatomically Constrained Tractography')
        brainIN=fullfile(OUTdir,'rT1_brain.nii.gz');
        file5tt=fullfile(OUTdir,'fsl5tt.nii.gz');
        
        sentence3=sprintf('5ttgen fsl -premasked %s %s',brainIN,file5tt);
      %  [~,fsl5tt_return]=system(sentence);
        
    % generate streamlines
  %  disp('4 - Generating Streamlines')
        fileStreamlines=fullfile(OUTdir,'10m_streamlines.tck');
        
        sentence4=sprintf('tckgen %s %s -act %s -crop_at_gmwmi -algorithm iFOD2 -seed_dynamic %s -select 10M',...
            fileFOD,fileStreamlines,file5tt,fileFOD);
     %   [~,tckgen_return]=system(sentence);
        
    % filter streamlines
 %   disp('5 - Running SIFT Filtering')
        fileFiltStreamlines=fullfile(OUTdir,'1m_sift_streamlines.tck');
        
        sentence5=sprintf('tcksift -force -act %s -term_number 1M %s %s %s',...
            file5tt,fileStreamlines,fileFOD,fileFiltStreamlines);
        %[~,sift_return]=system(sentence);
        
    % extract fronto-striatal connectivity
    disp('6 - Extracting FrontoStriatal Matrix')
        %register template to DWI space
        fileSeedImage='/datay2/chumin-F31/20190302_PET-DTI_final/seed_roi/frontostriatal.nii.gz';
        MNI2T1_linear=fullfile(T1dir,'registration/MNI2T1_linear.mat');
                if ~exist(MNI2T1_linear,'file')
                    dof6_2T1=fullfile(T1dir,'registration/MNI2T1_dof6.mat');
                    dof12_2T1=fullfile(T1dir,'registration/MNI2T1_dof12.mat');
                    sentence=sprintf('convert_xfm -omat %s -concat %s %s',...
                        MNI2T1_linear,dof12_2T1,dof6_2T1);
                    [~,result]=system(sentence);
                end
                % apply MNI2T1 transformation
                MNI2T1_warp=fullfile(T1dir,'registration/MNI2T1_warp.nii.gz');
                fileT1ref=fullfile(T1dir,'T1_fov_denoised.nii');
                [~,name,ext]=fileparts(fileSeedImage);
                T1seedImage=fullfile(OUTdir,sprintf('T1_%s',[name ext]));
                sentence=sprintf('applywarp -i %s -r %s -o %s -w %s --postmat=%s --interp=nn',...
                    fileSeedImage,fileT1ref,T1seedImage,MNI2T1_warp,MNI2T1_linear);
                [~,result]=system(sentence);
                % apply T1 to DWI transformation
                fileMat2FA=fullfile(DWIdir,'T1_2FAdof6.mat');
                filewarp=fullfile(DWIdir,'T12FA_warpfield.nii.gz');
                DWIseedImage=fullfile(OUTdir,sprintf('DWI_%s',[name ext]));
                fileRef=fullfile(DWIdir,'DTIfit/3_DWI_FA.nii.gz');
                fileT1brain=fullfile(T1dir,'T1_brain.nii.gz');
                DWI_T1brain=fullfile(OUTdir,'rT1_brain.nii.gz');
                fileRef1mm=fullfile(OUTdir,'/3_DWI_FA_1mm.nii.gz');
                sentence=sprintf('flirt -in %s -ref %s -out %s -applyisoxfm 1',...
                    fileRef,fileRef,fileRef1mm);
                [~,result]=system(sentence);
                if exist(filewarp,'file')
                    sentence=sprintf('applywarp -i %s -o %s -r %s -w %s --interp=nn',...
                        T1seedImage,DWIseedImage,fileRef,filewarp);
                    [~,result]=system(sentence);
                    sentence=sprintf('applywarp -i %s -o %s -r %s -w %s --interp=nn',...
                        fileT1brain,DWI_T1brain,fileRef1mm,filewarp);
                    [~,result]=system(sentence);
                else
                    sentence=sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                        T1seedImage,fileRef,fileMat2FA,DWIseedImage);
                    [~,result]=system(sentence);
                    sentence=sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour',...
                        fileT1brain,fileRef1mm,fileMat2FA,DWI_T1brain);
                    [~,result]=system(sentence);
                end

                % mask & dilate image within gray matter
                GM_Mask=fullfile(DWIdir,'rT1_GM_mask.nii.gz');
                fileMasked=fullfile(OUTdir,sprintf('DWI_GM_%s',[name ext]));
                sentence=sprintf('fslmaths %s -mas %s %s',DWIseedImage,GM_Mask,fileMasked);
                [~,result]=system(sentence);
                
            fileConnectome=fullfile(OUTdir,'1M_2mm_radial_frontostraital_connectome.csv');
                
            sentence6=sprintf('tck2connectome -assignment_radial_search 2 -symmetric -zero_diagonal -force %s %s %s',...
                fileFiltStreamlines,fileMasked,fileConnectome);
          %  [~,connectome_return]=system(sentence);
         
        fileout=fullfile(OUTdir,'run_mrtrix.sh');
        fID=fopen(fileout,'w');
        fprintf(fID,'#!/bin/bash\n%s\n%s\n%s\n%s\n%s\n%s\n',sentence1,sentence2,sentence3,sentence4,sentence5,sentence6);
        fclose(fID);
        system(sprintf('chmod +x %s',fileout))
%dwi2response tournier -voxels test_voxels.nii.gz -force -mask EDDY/b0_brain_mask.nii.gz -fslgrad EDDY/eddy_output.eddy_rotated_bvecs DTIfit/3_DWI.bval EDDY/eddy_output.nii.gz test_response.txt
%dwi2fod csd -fslgrad EDDY/eddy_output.eddy_rotated_bvecs DTIfit/3_DWI.bval -mask EDDY/b0_brain_mask.nii.gz EDDY/eddy_output.nii.gz test_response.txt test_fod.mif
%5ttgen fsl -premasked rT1_brain.nii.gz test_5tt.nii.gz
%tckgen test_fod.mif 10m.tck -act test_5tt.nii.gz -crop_at_gmwmi -algorithm iFOD2 -seed_dynamic test_fod.mif -select 10M 
%tcksift -force -act test_5tt.nii.gz -term_number 1M 10m.tck test_fod.mif 1m.tck
%tck2connectome -assignment_radial_search 2 -symmetric -zero_diagonal -force -keep_unassigned 1m.tck rT1_GM_parc_shen.nii.gz 1m_radial2mm_connectome_keptOrphans.csv 
    end
    clear sample   
end