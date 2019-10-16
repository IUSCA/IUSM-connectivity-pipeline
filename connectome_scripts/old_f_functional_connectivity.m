function [paths,flags,configs,parcs,params]=f_functional_connectivity(paths,flags,configs,parcs)
%                       F_FUNCTIONAL_CONNECTIVITY
% Functional MNI preprocessing, registration to anatony and preparation of
% connectivity matrices.
%                       
% Contributors:
%   Joaquin Goni, Purdue University
%   Joey Contreras, University Southern California
%   Enrico Amico, Purdue University
%   Mario Dzemidzic, Indiana University School of Medicine
%   Evgeny Chumin, Indiana University School of Medicine
%   Zikai Lin, Indiana University School of Medicine
%
%   Update:
%       03/26/2019 Integrated MELODIC and ICA-AROMA. -ZL 
%       04/30/2019 PCA not hardwired to 0 1 3 5 any more, can be defined by
%       user in batch_set_up.m

%% Convert dcm2nii
if flags.EPI.dcm2niix == 1
    disp('---------------------------------')
    disp('0. Dicom to NIFTI conversion')
    paths.EPI.dcm = fullfile(paths.EPI.dir,configs.name.dcmFolder);
    % Remove existing nifti image files
    epifile = '0_epi';
    fileNii = fullfile(paths.EPI.dir,sprintf('%s.nii',epifile));
    fileNiigz = fullfile(paths.EPI.dir,sprintf('%s.nii.gz',epifile));
    if exist(fileNii,'file') || exist(fileNiigz,'file')
        sentence = sprintf('rm -fr %s/%s.nii*',paths.EPI.dir,epifile);
        [~,result] = system(sentence); %#ok<*ASGLU>
    end
    % import dicoms
    fileLog= sprintf('%s/dcm2niix.log',paths.EPI.dir);
    sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.scripts,epifile,paths.EPI.dir,paths.EPI.dcm,fileLog);
    [~,result] = system(sentence);
    sentence=sprintf('gzip -f %s/0_epi.nii',paths.EPI.dir);
    [~,result] = system(sentence);
    
    if exist(fileNiigz,'file')
        disp('  - 0_epi.nii.gz file created')
        disp('---------------------------------')
    else
        fprintf('%s not created. Exiting...\n',fileNiigz)
        return
    end
else
    if exist(fullfile(paths.EPI.dir,'0_param_dcm_hdr.mat'),'file')
        load(fullfile(paths.EPI.dir,'0_param_dcm_hdr.mat'),'params');
    end
end



%% Read info from the headers of the dicom fMRI volumes
if flags.EPI.ReadHeaders==1
    paths.EPI.dcm = fullfile(paths.EPI.dir,configs.name.dcmFolder);
    disp('---------------------------------')
    disp('Dicom Header Information')
%-------------------------------------------------------------------------%
    % Identify dicoms
    dcmext = strcat('*.',configs.name.dcmFiles); % search for '*.dcm' or '*.IMA'
    dicom_files = dir(fullfile(paths.EPI.dcm,dcmext));
    if size(dicom_files,1) == 0
        warning('No dicom (.IMA or .dcm) images found. Skipping further analysis')
        return
    end
%-------------------------------------------------------------------------%
    if configs.EPI.UseJson == 1
        disp('Using json file provided by dcm2niix to extract header information...')
        % Read Json files from dcm2niix output directory
        dcm2niix_json_loc = fullfile(paths.EPI.dir, '0_epi.json');
        dcmHeaders_json = get_features_json(dcm2niix_json_loc,  true, true);
        
        params.EPI.TR = dcmHeaders_json.RepetitionTime; % TR
        fprintf(' -JSON: Repetition Time (TR): %d\n',dcmHeaders_json.RepetitionTime)
        
        params.EPI.TE = dcmHeaders_json.EchoTime; % TE
        fprintf(' -JSON: Echo Time (TE): %d\n',dcmHeaders_json.EchoTime)
        
        params.EPI.FlipAngle = dcmHeaders_json.FlipAngle; % Flip Angle
        fprintf(' -JSON: Flip Angle: %d\n',dcmHeaders_json.FlipAngle)
        
        params.EPI.EffectiveEchoSpacing = dcmHeaders_json.EffectiveEchoSpacing; % Effective Echo Spacing
        fprintf(' -JSON: Effective Echo Spacing: %g\n',dcmHeaders_json.EffectiveEchoSpacing)
        
        params.EPI.BandwidthPerPixelPhaseEncode = dcmHeaders_json.BandwidthPerPixelPhaseEncode; % Band width Per Pixel Phase Encode
        fprintf(' -JSON: Band width Per Pixel Phase Encode: %f\n',dcmHeaders_json.BandwidthPerPixelPhaseEncode)
        
        % Extract Slicing time
        params.EPI.slice_fractimes = dcmHeaders_json.SliceTiming;
        params.EPI.n_slice = size(dcmHeaders_json.SliceTiming, 1); % number of slice
        fprintf(' -DATA: Number of Slices: %d\n',params.EPI.n_slice)
        slice_fractimes_uniq = unique(params.EPI.slice_fractimes);
        if size(slice_fractimes_uniq,1) == size(params.EPI.slice_fractimes,1) % each slice acquired at a different time
            if issorted(params.EPI.slice_fractimes) % Sequential
                params.EPI.slice_ord = 1; 
                params.EPI.slice_rev = 0; % increasing from bottom to top (default)
            elseif issorted(params.EPI.slice_fractimes(params.EPI.n_slice:-1:1)) % Sequential
                params.EPI.slice_ord = 1; 
                params.EPI.slice_rev = 1; % increasing from top to bottom ("--down"")
            else
                % Interleaved?
                oddslice_fractimes = params.EPI.slice_fractimes(1:2:params.EPI.n_slice);
                n_oddslices = size(oddslice_fractimes,1);
                evenslice_fractimes = params.EPI.slice_fractimes(2:2:params.EPI.n_slice);
                n_evenslices = size(evenslice_fractimes,1);
                if max(oddslice_fractimes) < min(evenslice_fractimes) 
                    params.EPI.slice_ord = 0; 
                    if issorted(oddslice_fractimes) % 1,3,5...2,4,6...
                        params.EPI.slice_rev = 0; % bottom-->top
                    elseif issorted(oddslice_fractimes(n_oddslices:-1:1)) % 31,29,27,...32,30,28,...
                        params.EPI.slice_rev = 1; % top-->bottom
                    end
                elseif max(evenslice_fractimes) < min(oddslice_fractimes)
                   params.EPI.slice_ord = 0; 
                   if issorted(evenslice_fractimes) %2,4,6,...1,3,5,...
                       params.EPI.slice_rev = 0; % bottom-->top
                   elseif issorted(evenslice_fractimes(n_evenslices:-1:1)) %32,30,...31,29,... 
                       params.EPI.slice_rev = 1; % top--> bottom
                   end
                else
                    error('Slice time information is not consistent with interleaved acquisition')
                end
            end
        else            % Multiband (multiple slices acquired at a same time)
            params.EPI.slice_ord = 2;
        end        
    else
        % Extract Repetition time (TR)
        % Dicom header information --> Flag 0018,0080 "Repetition Time"
        [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0018,0080 ' paths.EPI.dcm '/' dicom_files(1).name]);
        params.EPI.TR = str2double(result(strfind(result,' '):end))/1000;
        fprintf(' -HEADER: Repetition Time (TR): %f\n',params.EPI.TR)
    %-------------------------------------------------------------------------%
        % Slice Time Acquisition
        dcm_file = fullfile(paths.EPI.dcm,dicom_files(1).name);
        sentence = sprintf('%s/dicom_hinfo -no_name -tag 0019,100a %s', paths.AFNI, dcm_file);
        [~,result] = system(sentence);
        %result
        n_slices = str2num(result); %#ok<*ST2NM>
        fprintf(' -HEADER: Number of Slices: %d\n',n_slices)

        id1 = 6;
        id2 = id1 + n_slices;
        sentence = sprintf('%s/dicom_hdr -slice_times %s |awk ''{split($0,a," "); for (i=%d;i<%d;i++) {print a[i]}}''', paths.AFNI, dcm_file,id1,id2);
        [~,result] = system(sentence);
        slice_times = str2double(strsplit(strtrim(result)))/1000;
        slice_fractimes = slice_times/params.EPI.TR;
        slicetimes_fname = fullfile(paths.EPI.dir,'slicetimes_frac.txt');
        fid = fopen(slicetimes_fname,'wt');
        fprintf(fid,'%f\n',slice_fractimes);
        fclose(fid);
        slice_fractimes_uniq = unique(slice_fractimes);
        if size(slice_fractimes_uniq,2) == size(slice_fractimes,2) % each slice acquired at a different time
            if issorted(slice_fractimes) % Sequential
                params.EPI.slice_ord = 1; 
                params.EPI.slice_rev = 0; % increasing from bottom to top (default)
            elseif issorted(slice_fractimes(n_slices:-1:1)) % Sequential
                params.EPI.slice_ord = 1; 
                params.EPI.slice_rev = 1; % increasing from top to bottom ("--down"")
            else
                % Interleaved?
                oddslice_fractimes = slice_fractimes(1:2:n_slices);
                n_oddslices = size(oddslice_fractimes,2);
                evenslice_fractimes = slice_fractimes(2:2:n_slices);
                n_evenslices = size(evenslice_fractimes,2);
                if max(oddslice_fractimes) < min(evenslice_fractimes) 
                    params.EPI.slice_ord = 0; 
                    if issorted(oddslice_fractimes) % 1,3,5...2,4,6...
                        params.EPI.slice_rev = 0; % bottom-->top
                    elseif issorted(oddslice_fractimes(n_oddslices:-1:1)) % 31,29,27,...32,30,28,...
                        params.EPI.slice_rev = 1; % top-->bottom
                    end
                elseif max(evenslice_fractimes) < min(oddslice_fractimes)
                   params.EPI.slice_ord = 0; 
                   if issorted(evenslice_fractimes) %2,4,6,...1,3,5,...
                       params.EPI.slice_rev = 0; % bottom-->top
                   elseif issorted(evenslice_fractimes(n_evenslices:-1:1)) %32,30,...31,29,... 
                       params.EPI.slice_rev = 1; % top--> bottom
                   end
                else
                    slicetime
                    error('Slice time information is not consistent with interleaved acquisition')
                end
            end
        else            % Multiband (multiple slices acquired at a same time)
            params.EPI.slice_ord = 2;
        end
        %% TE
        [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0018,0081 ' paths.EPI.dcm '/' dicom_files(1).name]);
        params.EPI.TE = str2double(result(strfind(result,' '):end));
        fprintf(' -HEADER: Echo Time (TE): %f\n',params.EPI.TE)
    end
%     %% esp
%     [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0043,102c ' paths.EPI.dcm '/' dicom_files(1).name]);
%     params.esp = str2double(result(strfind(result,' '):end));
%     fprintf('Header extracted esp: %f\n',params.esp)
%     %% asset
%     [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0043,1083 ' paths.EPI.dcm '/' dicom_files(1).name]);
%     tmp = result(strfind(result,' '):end);
%     params.asset = str2double(tmp(1:(strfind(tmp,'\')-1)));
%     fprintf('Header extracted asset: %f\n',params.asset)
    save(fullfile(paths.EPI.dir,'0_param_dcm_hdr.mat'),'params');
    fprintf('---------------------------------\n')
end

%% 0. Spin Echo Unwarping
if flags.EPI.SpinEchoUnwarp==1
    disp('-----------------------')
    disp('0. Field Map Correction')
    disp('-----------------------')
    
    % set up directory paths
    paths.EPI.SEFM = fullfile(paths.EPI.dir,configs.name.sefmFolder);
    paths.EPI.APdcm = fullfile(paths.EPI.SEFM,configs.name.APdcm);
    paths.EPI.PAdcm = fullfile(paths.EPI.SEFM,configs.name.PAdcm);
    paths.EPI.GREmagdcm = fullfile(paths.EPI.SEFM,configs.name.GREmagdcm);
    paths.EPI.GREphasedcm = fullfile(paths.EPI.SEFM,configs.name.GREphasedcm);
    
    if ~exist(paths.EPI.SEFM,'dir') 
        warning('%s',paths.EPI.SEFM,' does not exist. Field map correction must be skipped.')
        return
    else
        if exist(paths.EPI.APdcm,'dir') && exist(paths.EPI.PAdcm,'dir')
            fileNiiAP= 'AP';
            sentence=sprintf('rm -fr %s/%s.nii*',paths.EPI.SEFM,fileNiiAP);
            [~,result] = system(sentence); % remove any existing .nii images
            fileLog= sprintf('%s/dcm2niix_AP.log',paths.EPI.SEFM);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.spripts,fileNiiAP,paths.EPI.SEFM,paths.EPI.APdcm,fileLog);
            [~,result] = system(sentence); % import AP fieldmaps

            fileNiiPA= 'PA';
            sentence=sprintf('rm -fr %s/%s.nii*',paths.EPI.SEFM,fileNiiPA);
            [~,result] = system(sentence); % remove any existing .nii images
            fileLog= sprintf('%s/dcm2niix_PA.log',paths.EPI.SEFM);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.scripts,fileNiiPA,paths.EPI.SEFM,paths.EPI.PAdcm,fileLog);
            [~,result] = system(sentence); % import PA fieldmaps
            
            sentence=sprintf('gzip -f %s/AP.nii %s/PA.nii',paths.EPI.SEFM,paths.EPI.SEFM);
            [~,result] = system(sentence); % gzip fieldmap volumes
            
            %Concatenate the AP then PA into single 4D image
            fileOut = fullfile(paths.EPI.SEFM,'sefield.nii.gz');
            if exist(fileOut,'file')
                sentence=sprintf('rm -fr %s',fileOut);
                [~,result] = system(sentence); % if it already exists; remove it.
            end
            fileInAP = fullfile(paths.EPI.SEFM,'AP.nii.gz');
            fileInPA = fullfile(paths.EPI.SEFM,'PA.nii.gz');
            sentence = sprintf('%s/fslmerge -tr %s %s %s %f',paths.FSL,fileOut,fileInAP,fileInPA,params.EPI.TR);
            [~,result]=system(sentence);
            
            
            % Generate an acqparams text file based on number of field maps.
            configs.EPI.SEreadOutTime = get_readout(paths,configs.name.dcmFiles);
            fprintf("SEreadOutTime: %f\n",configs.EPI.SEreadOutTime);
            APstr=[0 -1 0 configs.EPI.SEreadOutTime];
            PAstr=[0 1 0 configs.EPI.SEreadOutTime];
            sentence = sprintf('%s/fslinfo %s |grep dim4|awk ''{split($0,a," "); {print a[2]}}''|head -n 1',paths.FSL, fileOut');
            [status,result] = system(sentence); % extract number of volumes from sefield.nii.gz
            
            if status == 0
                if rem(str2num(result),2) == 0
                    SEnumMaps = str2num(result)/2; % convert number of volumes to a number
                else
                    warning('sefile.nii.gz file must contain even number of volumes. Exiting...')
                    return
                end
            else
                SEnumMaps = configs.EPI.SEnumMaps; % trust the user input if SEnumMaps failed
            end
            
            for i=1:SEnumMaps
                acqparams(i,:)=APstr; %#ok<*AGROW>
            end
            for i=SEnumMaps+1:SEnumMaps*2
                acqparams(i,:)=PAstr;
            end
            fileParams = fullfile(paths.EPI.SEFM,'acqparams.txt');
            dlmwrite(fileParams,acqparams,'delimiter','\t')
            
            % Generate (topup) and apply (applytopup) spin echo field map
            % correction to 0_epi image.
            
            fileIn = fullfile(paths.EPI.SEFM,'sefield.nii.gz');
            if exist(fileParams,'file') && exist(fileIn,'file')
                fileOutName = fullfile(paths.EPI.SEFM,'topup_results');
                fileOutField = fullfile(paths.EPI.SEFM,'topup_field');
                fileOutUnwarped = fullfile(paths.EPI.SEFM,'topup_unwarped');
                if flags.EPI.RunTopup == 1
                    disp('Starting topup on sefield.nii.gz. This might take a while...')
                    sentence = sprintf('%s/topup --imain=%s --datain=%s --out=%s --fout=%s --iout=%s',...
                        paths.FSL,fileIn,fileParams,fileOutName,fileOutField,fileOutUnwarped);
                    [~,result]=system(sentence);
                end
                if exist(sprintf('%s.nii.gz',fileOutUnwarped),'file') ~= 2
                    warning('Topup output not created. Exiting...')
                    return
                end
                
                fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
                if exist(fileIn,'file')
                    fileOut = fullfile(paths.EPI.SEFM,'0_epi_unwarped.nii.gz');
                    sentence = sprintf('%s/applytopup --imain=%s --datain=%s --inindex=1 --topup=%s --out=%s --method=jac',...
                        paths.FSL,fileIn,fileParams,fileOutName,fileOut);
                    disp('Starting applytopup on 0_epi.nii.gz.')
                    [~,result]=system(sentence);
                else
                    warning('0_epi.nii.gz not found. Exiting...')
                    return
                end
                if exist(fullfile(paths.EPI.SEFM,'0_epi_unwarped.nii.gz'),'file')
                    sentence = sprintf('mv %s %s',fullfile(paths.EPI.SEFM,'0_epi_unwarped.nii.gz'),fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz'));
                    [status,result] = system(sentence);
                    if status == 0
                        disp('EPI volume unwarping completed.')
                        fprintf('---------------------------------\n')
                    end
                end
            else
                warning('UNWARP/sefield.nii.gz or acqparams.txt are missing. topup not started')
                return
            end
        elseif exist(paths.EPI.GREmagdcm,'dir') && exist(paths.EPI.GREphasedcm,'dir')
            % Identify dicoms
            dcmext = strcat('*.',configs.name.dcmFiles); % search for '*.dcm' or '*.IMA'
            dicom_files = dir(fullfile(paths.EPI.GREmagdcm,dcmext));
            if size(dicom_files,1) < 2
                warning('No dicom (.IMA or .dcm) images found. Skipping further analysis')
                return
            else
                % Extract TE1 and TE2 from the first image of Gradient Echo Magnitude Series
                % fsval image descrip would do the same but truncates TEs to a single digit!
                [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0018,0081 ' paths.EPI.GREmagdcm '/' dicom_files(1).name]);
                aux = strfind(result,' ');
                TE1 = str2num(result(aux+1:end));
                [~,result] = system([paths.AFNI '/' 'dicom_hinfo -tag 0018,0081 ' paths.EPI.GREmagdcm '/' dicom_files(2).name]);
                aux = strfind(result,' ');
                TE2 = str2num(result(aux+1:end));
                fprintf('Header extracted TE1 = %5.2f ms, TE2 = %5.2f ms\n',TE1,TE2)
                DeltaTE = TE2 - TE1;   
            end
            fileNm1 = 'mag';
            % remove any existing images
            fileMag1 = fullfile(paths.EPI.SEFM,'mag.nii');
            if exist(fileMag1,'file'); delete(fileMag1); end
            fileMag2 = fullfile(paths.EPI.SEFM,'_e2mag.nii');
            if exist(fileMag2,'file'); delete(fileMag2); end 
            % dicom import
            fileLog= sprintf('%s/dcm2niix.log',paths.EPI.GREmagdcm);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.scripts,fileNm1,paths.EPI.SEFM,paths.EPI.GREmagdcm,fileLog);
            [~,result] = system(sentence);
            
            fileNm1 = '_e2'; % dcm2niix automatically prepends for 2nd echo images
            fileNm2 = 'gre_fieldmap_phasediff';
            filePhaseMap1 = fullfile(paths.EPI.SEFM,strcat(fileNm1,fileNm2,'.nii'));
            filePhaseMap2 = fullfile(paths.EPI.SEFM,strcat(fileNm2,'.nii'));
            % remove any existing file
            if exist(filePhaseMap1,'file'); delete(filePhaseMap1); end
            % dicom import
            fileLog=sprintf('%s/dcm2niix.log',paths.EPI.GREphasedcm);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',...
                paths.scripts,fileNm2,paths.EPI.SEFM,paths.EPI.GREphasedcm,fileLog);
            [~,result] = system(sentence);
            % Copy and gzip the nifti images
            if exist(filePhaseMap1,'file')
                sentence = sprintf('cp %s %s',filePhaseMap1,filePhaseMap2);
                [~,result] = system(sentence);
                sentence = sprintf('gzip -f %s',filePhaseMap2);
                [~,result] = system(sentence);
                filePhaseMap = sprintf('%s.gz',filePhaseMap2);
            else
                warning('%s not found! Exiting...',filePhaseMap1)
                return
            end
            
            if exist(fileMag1,'file') && exist(fileMag2,'file')
                fileMagAvg = fullfile(paths.EPI.SEFM,'gre_fieldmap_magAVG');
                sentence = sprintf('%s/fslmaths %s -add %s -div 2 %s',paths.FSL,...
                    fileMag1,fileMag2,fileMagAvg);
                [~,result] = system(sentence);
                
                fileIn = sprintf('%s.nii.gz',fileMagAvg);
                fileMagBrain = fullfile(paths.EPI.SEFM,'gre_fieldmap_magAVG_brain');
                sentence = sprintf('%s/bet %s %s -f %.4f -g %.4f -m',paths.FSL,...
                    fileIn,fileMagBrain,configs.EPI.GREbetf,configs.EPI.GREbetg);
                [~,result] = system(sentence);
                
                % Prepare phase map
                fileFMap = fullfile(paths.EPI.SEFM,strcat(fileNm2,'_prepared'));
                sentence = sprintf('%s/fsl_prepare_fieldmap SIEMENS %s %s %s %.4f',...
                    paths.FSL,filePhaseMap,fileMagBrain,fileFMap, DeltaTE);
                [~,result] = system(sentence);
 
                % Now run fugue (-s 3 : apply Gaussian smoothing of sigma = 3 mm
                fileFMapIn = sprintf('%s.nii.gz',fileFMap);

                if configs.EPI.GREdespike == 1
                    fileFMapRads = fullfile(paths.EPI.SEFM,sprintf('fm_rads_brain_sm%d_m_ds',configs.EPI.GREsmooth));
                    sentence = sprintf('%s/fugue --loadfmap=%s -s %.2f -m --despike --savefmap=%s',...
                        paths.FSL,fileFMapIn,configs.EPI.GREsmooth,fileFMapRads);
                else
                    fileFMapRads = fullfile(paths.EPI.SEFM,sprintf('fm_rads_brain_sm%d_m',configs.EPI.GREsmooth));
                    sentence = sprintf('%s/fugue --loadfmap=%s -s %.2f -m --savefmap=%s',...
                        paths.FSL,fileFMapIn,configs.EPI.GREsmooth,fileFMapRads);
                end
                [~,result] = system(sentence);
                
                fileFMapRadsOut = sprintf('%s.nii.gz',fileFMapRads);
                if exist(fileFMapRadsOut,'file')
                    disp('Fugue successfully created fm_rads_brain field map.')
                else
                    warning('%s not created: fugue failed. Exiting...',fileFMapRadOut)
                    return
                end
% Now upwarp 0_epi.nii.gz
% Echo Spacing = 1 / [(0019,1028) * (0051,100B component #1)]
% where 0019,1028 is 40.584
% and (0051,100B) is 80 per [80*80]
% SoL denominator is 3246.72 Hz
% Echo Spacing = 0.000308
% Echo Spacing Siemens PDF="0.000616"/GRAPPA="2" = 0.000308 
% TIME="0.000308" sec
% fugue -i f$STUDY$EXAM-$index2.nii --dwell=$DWELLTIME --loadfmap=$FMAPDIR/$fmap --unwarpdir=y- -u uf$STUDY$EXAM-$index2

                fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
                if exist(fileIn,'file')
                    fileOut = fullfile(paths.EPI.dir,'0_epi_unwarped');
                    sentence = sprintf('%s/fugue -i %s --dwell=%.6f --loadfmap=%s --unwarpdir=y- -u %s',...
                        paths.FSL,fileIn,configs.EPI.EPIdwell,fileFMapRadsOut,fileOut);
                    disp('Applying fugue to unwarp 0_epi.nii.gz.')
                    [~,result]=system(sentence);
                    if exist(sprintf('%s.nii.gz',fileOut),'file')
                        disp('0_epi_unwarped.nii.gz successfully created.')
                    else
                        warning('fugue unwarping failed. Exiting...')
                    end
                else
                    warning('0_epi.nii.gz not found. Exiting...')
                    return
                end
                
            else
                warning('Field Map images not created. Exiting...')
                return
            end
            
        else
            warning('UNWARP DICOMS folders do exist. Field Map correction failed.')
            return
        end
    end
end

%% 1. Slice timing correction
if flags.EPI.SliceTimingCorr==1
    
    disp('------------------------------------')
    disp('1. Slice Time Acquisition Correction')
    disp('------------------------------------')
    
    if ~exist(fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
    else
        if flags.EPI.UseUnwarped == 1
            fileIn = fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz');
            disp(' -Processing: 0_epi_unwarped.nii.gz')
        else
            fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
            disp(' -Processing: 0_epi.nii.gz')
        end
    end

    if ~exist(fileIn,'file') % check that selected EPI exists.
        warning('File %s',fileIn,'does not exist. Skipping further analysis')
        return
    end
    sentence = sprintf('%s/fslreorient2std %s %s',paths.FSL,fileIn,fileIn);
    [~,result] = system(sentence);
    
    fileRef = fullfile(paths.EPI.dir,'slicetimes_frac.txt');
    fileOut = fullfile(paths.EPI.dir,'1_epi.nii.gz');
    if params.slice_ord == 2 || configs.EPI.UseTcustom == 1 % Use custom interleave timing file
        sentence = sprintf('%s/slicetimer -i %s -o %s -r %0.4f --tcustom=%s',...
            paths.FSL,fileIn,fileOut,params.TR,fileRef);
    elseif params.slice_ord == 1 && configs.EPI.UseTcustom ~= 1 % Sequential acquisition
        if params.slice_rev == 0
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f';
        elseif params.slice_rev == 1
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --down';
        end
        sentence = sprintf(st_command,paths.FSL,fileIn,fileOut,params.TR);
    elseif params.slice_ord == 0 && configs.EPI.UseTcustom ~= 1 % Interleaved acquisition
        if params.slice_rev == 0
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --odd';
        elseif params.slice_rev == 1
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --odd --down';
        end
        sentence = sprintf(st_command,paths.FSL,fileIn,fileOut,params.TR);
    end
    [~,result] = system(sentence);
end

%% 2. Motion correction
if flags.EPI.MotionCorr==1
    disp('--------------------')
    disp('2. Motion Correction')
    disp('--------------------')
    
    if ~exist(fullfile(paths.EPI.dir,'1_epi.nii.gz'),'file')
        disp(' -No slice time corrected 1_epi output found;')
        disp('  Defaulting to 0_epi data.')
        if ~exist(fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz'),'file')
            disp(' -Unwarped 0_epi volume does not exist')
            fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
            if exist(fileIn,'file')
                disp(' -Will use 0_epi from dicom conversion.')
            else
                warning(' -No 0_epi inputs found exiting...')
                return
            end
        else
            if flags.EPI.UseUnwarped == 1
                fileIn = fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz');
                disp(' -Will use 0_epi_unwarped.nii.gz as set by UseUnwarped flag =1')
            else
                fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
                disp(' -Will use 0_epi.nii.gz as set by UseUnwarped flag =0')
            end
        end
    else
        fileIn = fullfile(paths.EPI.dir,'1_epi.nii.gz');
        disp(' -Will use the slice time corrected 1_epi.nii.gz as input')
    end
      
        fileOut = fullfile(paths.EPI.dir,'2_epi');
        sentence=sprintf('%s/fslval %s dim4',paths.FSL,fileIn);
        [~,result] = system(sentence);
        params.EPI.nvols = str2double(result);
        fprintf('Number of volumes in 1_epi_brain: %d\n',params.EPI.nvols)

        sentence = sprintf('%s/mcflirt -in %s -out %s -plots -meanvol',...
            paths.FSL,fileIn,fileOut);
        [~,result] = system(sentence);
        
        sentence = sprintf('mv %s/2_epi.par %s/motion.txt',paths.EPI.dir,paths.EPI.dir);
        [~,result] = system(sentence);
end

%% 3 bet fmri and T1 registration
if flags.EPI.RegT1==1
    disp('-------------------------------')
    disp('3. BEt fMRI and T1 Registration')
    disp('-------------------------------')
    
    if exist(fullfile(paths.EPI.dir,'2_epi.nii.gz'),'file') == 0
        warning('File %s',fullfile(paths.EPI.dir,'2_epi.nii.gz'),'does not exist. Skipping further analysis')
        return
    end
%-------------------------------------------------------------------------%    
    % compute the mean-vol of epi along 4th dimension (time)
    fileIn = fullfile(paths.EPI.dir,'2_epi.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'2_epi_meanvol.nii.gz');
    sentence = sprintf('%s/fslmaths %s -Tmean %s',paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    sentence = sprintf('%s/bet %s %s -f %.4f -n -m -R',paths.FSL,fileOut,fileOut,configs.EPI.epibetF);
    [~,result] = system(sentence);
    
    fileIn = fullfile(paths.EPI.dir,'2_epi.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'3_epi.nii.gz');
    fileMas = fullfile(paths.EPI.dir,'2_epi_meanvol_mask.nii.gz');
    sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileIn,fileMas,fileOut);
    [~,result] = system(sentence);
%-------------------------------------------------------------------------%    
    % rigid body registration (dof 6) registration of T1 to fMRI
    fileIn = fullfile(paths.T1.dir,'T1_brain.nii.gz');
    fileRef = fullfile(paths.EPI.dir,'2_epi_meanvol.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_brain_dof6');
    fileOmat = fullfile(paths.EPI.dir,'T1_2_epi_dof6.mat');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -omat %s -cost normmi -dof 6 -interp spline',...
        paths.FSL,fileIn,fileRef,fileOut,fileOmat);
    [~,result] = system(sentence);
    % generate an inverse transform fMRI to T1
    fileOmatInv = fullfile(paths.EPI.dir,'epi_2_T1_dof6.mat');
    sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',...
        paths.FSL,fileOmatInv,fileOmat);
    [~,result] = system(sentence);
%-------------------------------------------------------------------------%    
    % apply transformation to T1 WM mask.
    fileIn = fullfile(paths.T1.dir,'T1_WM_mask.nii.gz');
    fileRef = fullfile(paths.EPI.dir,'2_epi_meanvol.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_WM_mask');
    fileInit = fullfile(paths.EPI.dir,'T1_2_epi_dof6.mat');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
        paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
%-------------------------------------------------------------------------%    
    % bbr registration of fMRI to rT1_dof6 based on WMseg
    fileRef = fullfile(paths.EPI.dir,'rT1_brain_dof6.nii.gz');
    fileIn = fullfile(paths.EPI.dir,'2_epi_meanvol.nii.gz');
    fileOmat = fullfile(paths.EPI.dir,'epi_2_T1_bbr.mat');
    fileWMseg = fullfile(paths.EPI.dir,'rT1_WM_mask');
    sentence = sprintf('%s/flirt -in %s -ref %s -omat %s -wmseg %s -cost bbr',...
        paths.FSL,fileIn,fileRef,fileOmat,fileWMseg);
    [~,result] = system(sentence);
    % Generate inverse matrix of bbr (T1_2_epi)
    fileMat = fullfile(paths.EPI.dir,'epi_2_T1_bbr.mat');
    fileMatInv = fullfile(paths.EPI.dir,'T1_2_epi_bbr.mat');
    sentence = sprintf('%s/convert_xfm -omat %s -inverse %s',...
        paths.FSL,fileMatInv,fileMat);
    [~,result] = system(sentence);
    % Join the T1_2_epi dof6 and bbr matrices    
    fileMat1 = fullfile(paths.EPI.dir,'T1_2_epi_dof6.mat');
    fileMat2 = fullfile(paths.EPI.dir,'T1_2_epi_bbr.mat');
    fileMatJoint = fullfile(paths.EPI.dir,'T1_2_epi_dof6_bbr.mat');
    sentence = sprintf('%s/convert_xfm -omat %s -concat %s %s',...
        paths.FSL,fileMatJoint,fileMat2,fileMat1);
    [~,result] = system(sentence);
    % Join the epi_2_T1 dof and bbr matrices
    fileMat1 = fullfile(paths.EPI.dir,'epi_2_T1_bbr.mat');
    fileMat2 = fullfile(paths.EPI.dir,'epi_2_T1_dof6.mat');
    fileMatJoint = fullfile(paths.EPI.dir,'epi_2_T1_bbr_dof6.mat');
    sentence = sprintf('%s/convert_xfm -omat %s -concat %s %s',...
        paths.FSL,fileMatJoint,fileMat2,fileMat1);
    [~,result] = system(sentence);
%-------------------------------------------------------------------------%
end

%% Apply transformation to tissue and parcellation images.
if flags.EPI.RegOthers==1
      if exist(fullfile(paths.EPI.dir,'T1_2_epi_dof6_bbr.mat'),'file') == 0
        warning('File %s',fullfile(paths.EPI.dir,'T1_2_epi_dof6_bbr.mat'),'does not exist. Skipping further analysis')
        return
      end
%-------------------------------------------------------------------------%    
    % brain mask
    fileIn = fullfile(paths.T1.dir,'T1_brain_mask_filled.nii.gz');
    fileRef = fullfile(paths.EPI.dir,'2_epi_meanvol_mask.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_brain_mask');
    fileInit = fullfile(paths.EPI.dir,'T1_2_epi_dof6_bbr.mat');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
    % WM mask
    fileIn = fullfile(paths.T1.dir,'T1_WM_mask.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_WM_mask');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence); 
    % Eroded WM_mask
    fileIn = fullfile(paths.T1.dir,'T1_WM_mask_eroded.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_WM_mask_eroded');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
    % CSF mask
    fileIn = fullfile(paths.T1.dir,'T1_CSF_mask.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_CSF_mask');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
    % Eroded CSF mask
    fileIn = fullfile(paths.T1.dir,'T1_CSF_mask_eroded.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_CSF_mask_eroded');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
    % CSF ventricle mask
    fileIn = fullfile(paths.T1.dir,'T1_CSFvent_mask_eroded.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_CSFvent_mask_eroded');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour -nosearch',...
    paths.FSL,fileIn,fileRef,fileOut,fileInit);
    [~,result] = system(sentence);
%-------------------------------------------------------------------------%
    %Probabilistic GM reg to fMRI space
    fileIn = fullfile(paths.T1.dir,'T1_GM_mask.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_GM_mask_prob');
    sentence = sprintf('%s/flirt -applyxfm -init %s -interp spline -in %s -ref %s -out %s -nosearch',...
    paths.FSL,fileInit,fileIn,fileRef,fileOut);
    [~,result] = system(sentence);
    % binarize GM probability map
    fileIn = fullfile(paths.EPI.dir,'rT1_GM_mask_prob.nii.gz');
    fileOut = fullfile(paths.EPI.dir,'rT1_GM_mask');
    sentence = sprintf('%s/fslmaths %s -thr %0.4f -bin %s',...
    paths.FSL,fileIn,configs.EPI.GMprobthr,fileOut);
    [~,result] = system(sentence);
%%
% Apllying T1 to EPI transformations to parcellations
    for k=1:length(parcs.pdir)
        if parcs.pnodal(k).true == 1 % treat as a parcelattion that will serve as noded for connectivity
            %transformation from T1 to epi space
            fileIn = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'_dil.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_parc_',parcs.plabel(k).name,'.nii.gz'));
            sentence = sprintf('%s/flirt -applyxfm -init %s -interp nearestneighbour -in %s -ref %s -out %s -nosearch',...
                paths.FSL,fileInit,fileIn,fileRef,fileOut);
            [~,result] = system(sentence);
            % masking Shen with GM
            fileIn = fullfile(paths.EPI.dir,strcat('rT1_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileMul = fullfile(paths.EPI.dir,'rT1_GM_mask.nii.gz');
            sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileIn,fileMul,fileOut);
            [~,result] = system(sentence);
            % removal of small clusters within ROIs
            fileIn = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'_clean.nii.gz'));
            f_get_largest_clusters_only(fileIn,fileOut,configs.EPI.minVoxelsClust);
        elseif parcs.pnodal(k).true == 0 % treat as an organizational parcellation to group nodes
            % Added by MDZ; 10/06/2015
            % transformation from T1 to epi space
            fileIn = fullfile(paths.T1.dir,strcat('T1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            sentence = sprintf('%s/flirt -applyxfm -init %s -interp nearestneighbour -in %s -ref %s -out %s -nosearch',...
                paths.FSL,fileInit,fileIn,fileRef,fileOut);
            [~,result] = system(sentence);
            % dilate Yeo7 parcellation
            fileIn = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            sentence = sprintf('%s/fslmaths %s -dilD %s',paths.FSL,fileIn,fileOut);
            [~,result] = system(sentence);
            % masking Yeo7 with GM
            fileIn = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileMul = fullfile(paths.EPI.dir,'rT1_GM_mask.nii.gz');
            sentence = sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileIn,fileMul,fileOut);
            [~,result] = system(sentence);
            % removal of small clusters within ROIs
            fileIn = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'.nii.gz'));
            fileOut = fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'_clean.nii.gz'));
            f_get_largest_clusters_only(fileIn,fileOut,configs.EPI.minVoxelsClust);
        else
            warning('parcs.pnodal.true is not specified for %s parcellation. Transformation to EPI not done',parcs.plabel(k).name)
        end
    end
end

%% 4. NORMALIZATION TO 4D MEAN OF 1000
if flags.EPI.IntNorm4D == 1
    disp('--------------------------------')
    disp('4. Normalization to 4D mean 1000')
    disp('--------------------------------')
    
    if exist(fullfile(paths.EPI.dir,'3_epi.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'3_epi.nii.gz');
        fileOut = fullfile(paths.EPI.dir,'4_epi.nii.gz');
        sentence = sprintf('%s/fslmaths %s -ing 1000 %s',paths.FSL,fileIn,fileOut);
        [~,result] = system(sentence);
    else
        warning('File %s',fullfile(paths.EPI.dir,'3_epi.nii.gz'),'does not exist. Skipping further analysis')
        return
    end
end

%% 5. ICA-AROMA - Denoising
if flags.EPI.AROMA == 1
    disp('-----------------------')
    disp('5. ICA-AROMA: Denoising')
    disp('-----------------------')
    
    if exist(fullfile(paths.EPI.dir,'4_epi.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'4_epi.nii.gz');
    else
        warning('File %s',fullfile(paths.EPI.dir,'4_epi.nii.gz'),'does not exist. Terminating preprocessing.')
        return
    end
    
    paths.EPI.AROMA = fullfile(paths.EPI.dir,'AROMA');
    if ~exist(paths.EPI.AROMA,'dir')
        mkdir(paths.EPI.AROMA)
    end
    paths.EPI.AROMAreg = fullfile(paths.EPI.AROMA,'registration');
    if ~exist(paths.EPI.AROMAreg,'dir')
        mkdir(paths.EPI.AROMAreg)
    end
    
    disp(' - Generating Inputs')
    % EPI -> T1 linear transformation
    disp('   - EPI to T1 linear transform')
    fileMat = fullfile(paths.EPI.dir,'epi_2_T1_bbr_dof6.mat');
    if ~exist(fileMat,'file')
        warning('Linear EPI->T1 tanformation not found. Please run the RegT1 flag. Exiting...')
        return
    end
    % T1 -> MNI152_2mm nonlinear transformation
    disp('   - T1 to MNI 2mm nonlinear transform')
    fileT1 = fullfile(paths.T1.dir,'T1_brain.nii');
    fileMNI2mm = fullfile(paths.MNIparcs,'MNI_templates/MNI152_T1_2mm_brain.nii.gz');
    filedof12mat = fullfile(paths.EPI.AROMAreg,'T1_2_MNI2mm_dof12.mat');
    filedof12img = fullfile(paths.EPI.AROMAreg,'rT1_dof12_2mm.nii.gz');
    if ~exist(filedof12mat,'file')
        sentence = sprintf('%s/flirt -in %s -ref %s -omat %s -dof 12 -cost mutualinfo -interp spline -out %s',...
            paths.FSL,fileT1,fileMNI2mm,filedof12mat,filedof12img);
        [~,result]=system(sentence);
    else
        disp('     - Using existing T1->MNI_2mm linear transformation')
    end
    fileWarpImg = fullfile(paths.EPI.AROMAreg,'rT1_warped_2mm.nii.gz');
    fileWarpField = fullfile(paths.EPI.AROMAreg,'T1_2_MNI2mm_warpfield.nii.gz');
    if ~exist(fileWarpField,'file')
        sentence = sprintf('%s/fnirt --in=%s --ref=%s --aff=%s --iout=%s --cout=%s',...
            paths.FSL,fileT1,fileMNI2mm,filedof12mat,fileWarpImg,fileWarpField);
        [~,result]=system(sentence);
    else
        disp('     - Using existing T1->MNI_2mm nonlinear transformation')
    end
    % 6mm FWHM EPI data smoothing
    disp('   - Smoothing EPI data by 6mm FWHM')
    fileEPI = fullfile(paths.EPI.dir,'4_epi.nii.gz');
    fileSmooth = fullfile(paths.EPI.AROMA,'s6_4_epi.nii.gz');
    if ~exist(fileSmooth,'file')
        sentence = sprintf('%s/fslmaths %s -kernel gauss 2.547965400864057 -fmean %s',...
            paths.FSL,fileEPI,fileSmooth);
        [~,result]=system(sentence);
    else
        disp('     - Using existing smoothed s6_4_epi data')
    end
    % mcFLIRT realignment parameters
    disp('   - mcFLIRT realignment parameters')
    fileMovePar = fullfile(paths.EPI.dir,'motion.txt');
    if ~exist(fileMovePar,'file')
        warning('Movement parameters from mcFLIRT not found. Please run the MotionCorr flag. Exiting...')
        return
    end
    disp(' - Starting ICA-AROMA')
    paths.EPI.AROMAout = fullfile(paths.EPI.AROMA,'AROMA-output');
    sentence = sprintf('%s %s/ICA-AROMA/ICA_AROMA.py -in %s -out %s -mc %s -affmat %s -warp %s -overwrite',...
        paths.python,paths.scripts,fileSmooth,paths.EPI.AROMAout,fileMovePar,fileMat,fileWarpField);
    [~,result]=system(sentence);
    disp(result)
    if ~exist(fullfile(paths.EPI.AROMAout,'denoised_func_data_nonaggr.nii.gz'),'file')
        warning(' -AROMA output file not found!')
        fprintf(' -Exiting... Check AROMA directory, it could be melodic ICA\n')
        fprintf('  found too many components and AROMA did not filter porperly\n')
        fprintf('  If that is the case fslfilt can be ran manually.\n')
    else    
        
        disp('   - Done.')
    end
end

%% 6. DEMEAN AND DETREND
if flags.EPI.DemeanDetrend == 1
    disp('---------------------')
    disp('6. Demean and Detrend')
    disp('---------------------')

    if exist(fullfile(paths.EPI.dir,'AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz');
        disp(' -Working on AROMA output volume')
    else
        disp(' -No AROMA output found.')
        if exist(fullfile(paths.EPI.dir,'4_epi.nii.gz'),'file')
            disp('   -Working on 4_epi volume')
            fileIn = fullfile(paths.EPI.dir,'4_epi.nii.gz');
        else
            warning('No AROMA or 4_epi volume exists. Exiting...')
            return
        end
    end      
    
    % read data
    resting = MRIread(fileIn);
    [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);

    % read brain mask
    volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask.nii.gz'));
    volRef = MRIread(fullfile(paths.EPI.dir,'2_epi_meanvol_mask.nii.gz'));
    volBrain.vol = (volBrain.vol>0) & (volRef.vol~=0);
    fileOut = fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz');
    MRIwrite(volBrain,fileOut,'double');
    % fill holes in the brain mask, without changing FOV
    sentence = sprintf('%s/fslmaths %s -fillh %s',paths.FSL,fileOut,fileOut);
    [~,result] = system(sentence);

    %vFit = (1:numTimePoints)';
    for i=1:sizeX
        for j=1:sizeY
            for k=1:sizeZ
                if volBrain.vol(i,j,k)>0
                    TSvoxel = reshape(resting.vol(i,j,k,:),[numTimePoints,1]);
                    % demean & detrend
                    TSvoxel_detrended = detrend(TSvoxel-mean(TSvoxel),'linear');
                    resting.vol(i,j,k,:) = TSvoxel_detrended;
                end
            end
        end
        if (mod(i,25)==0)
            disp(i/sizeX)
        end
    end
    fileOut2 = fullfile(paths.EPI.dir,'6_epi.nii.gz');
    MRIwrite(resting,fileOut2,'double');
    sentence=sprintf('%s/fslmaths %s -mas %s %s',paths.FSL,fileOut2,fileOut,fileOut2);
    [~,result]=system(sentence);
end

%% 7. MOTION AND OUTLIER REGRESSORS
if flags.EPI.MotionRegressors == 1
    fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
    if exist(fileIn,'file')
        sentence = sprintf('%s/fslval %s pixdim4',paths.FSL,fileIn);
        [~,result] = system(sentence);
        params.TR = str2double(result);
        dropTRs = ceil(configs.EPI.scrubtime/params.TR);
        sentence = sprintf('%s/fslval %s dim4',paths.FSL,fileIn);
        [~,result] = system(sentence);
        params.nvols = str2double(result);
        TR.init = dropTRs+1;
        TR.end = params.nvols -dropTRs;
    else
        warning('File %s',fileIn,'does not exist. Skipping further analysis')
        return
    end
    
    if flags.EPI.GS ==1 
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
        GSreg_name = 'GSreg\_yes';
        if ~exist(paths.EPI.epiGS,'dir')
            mkdir(paths.EPI.epiGS);
        end
    elseif flags.EPI.GS ==0
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
        GSreg_name = 'GSreg\_no';
        if ~exist(paths.EPI.epiGS,'dir')
            mkdir(paths.EPI.epiGS);
        end
    else
        warning('flags.EPI.GS not specified. Exiting...')
    end

    disp('--------------------------------')
    disp('7. Motion and Outlier Regressors')
    disp('--------------------------------')

 if exist(fullfile(paths.EPI.dir,'6_epi.nii.gz'),'file')
    resting = MRIread(fullfile(paths.EPI.dir,'6_epi.nii.gz'));
    [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);
    if numTimePoints ~= params.nvols
        warning('Number of time points is inconsistent. Exiting...')
        return
    end
 else
     warning('File %s',fullfile(paths.EPI.dir,'6_epi.nii.gz'),'does not exist. Skipping further analysis')
        return
 end
%-------------------------------------------------------------------------%
    % load 6 motion regressors
    load(fullfile(paths.EPI.dir,'motion.txt')); %#ok<*LOAD>
    % derivatives of 6 motion regressors
    motion_deriv = nan(size(motion)); %#ok<*NODEF>
    for i=1:size(motion,2)
        m = motion(:,i)'; %#ok<*IDISVAR>
        m_deriv = [0,diff(m)];
        motion_deriv(:,i)=m_deriv';
    end
%-------------------------------------------------------------------------%
    fileMask = fullfile(paths.EPI.dir,'1_epi_brain_mask.nii.gz');
    % Frame Displacement regressor
    disp('Computing fd regressors')
    % example: fsl_motion_outliers -i 4_drafREST1.nii.gz -o motionOutputFD.txt -s motionMetricsFD.txt -p motionMetricsFD.png --fd
    metric = 'fd';
    fileIn = fullfile(paths.EPI.dir,'1_epi.nii.gz'); %only after slice-timing correction
    fileOut = fullfile(paths.EPI.epiGS,sprintf('motionRegressor_%s.txt',metric));
    fileMetrics = fullfile(paths.EPI.epiGS,sprintf('motionMetric_%s.txt',metric));
    if exist(fileOut,'file')
        delete(fileOut);
    end
     if exist(fileMetrics,'file')
        delete(fileMetrics);
    end
    filePlot = fullfile(paths.EPI.epiGS,sprintf('motionPlot_%s.png',metric));
    sentence = sprintf('%s/fsl_motion_outliers -i %s -o %s -s %s -p %s --%s --thresh=%0.2f -m %s',...
        paths.FSL,fileIn,fileOut,fileMetrics,filePlot,metric,configs.EPI.FDth,fileMask);
    [~,result] = system(sentence);
   
    if exist(fileMetrics,'file') 
        load(fileMetrics); % motionMetric_fd
        if ~isempty(configs.EPI.FDth) && ~isnan(configs.EPI.FDth)
            scrubbing_fd = motionMetric_fd < configs.EPI.FDth;
%             if nnz(~scrubbing_fd) >= 64
%                 fprintf('-- identified %d FD outliers\n',nnz(~scrubbing_fd))
%                 error('Too many FD outliers being scrubbed.')
%             end
            fprintf('-- identified %d FD outliers\n',nnz(~scrubbing_fd))
        else
            error('configs.EPI.FDth is undefined!')
        end
    else
        error('file %s not found!',fileMetrics)
    end
%-------------------------------------------------------------------------%
    % Standardized Dvars
    disp('Computing dvars regressors')
    metric = 'standardized_dvars';
    fileIn = fullfile(paths.EPI.dir,'1_epi.nii.gz');
    fileMetrics = fullfile(paths.EPI.epiGS, sprintf('motion_metric_%s.txt',metric));
    if exist(fileMetrics,'file')
        delete(fileMetrics);
    end
    
    % Make sure the DVARS.sh is executable
    sentence = sprintf('chmod +x %s', fullfile(paths.scripts, 'DVARS.sh'));
    system(sentence);
    
    % Now generate fileMetrics
    sentence = sprintf('%s %s %s',fullfile(paths.scripts,'DVARS.sh'),fileIn,fileMetrics);
    [~,result]=system(sentence);
    
    if exist(fileMetrics,'file')
        stDVARS = dlmread(fileMetrics,'\t',[0 0 numTimePoints-2 0]);
        stDVARS(2:end+1)=stDVARS; stDVARS(1) = 0; 
            if isempty(configs.EPI.DVARSth) || isnan(configs.EPI.DVARSth)
                error('configs.EPI.DVARSth is undefined!') 
            else
                scrubbing_dvars = stDVARS < configs.EPI.DVARSth;
            end
        fprintf('-- identified %d DVARS outliers\n',nnz(~scrubbing_dvars))       
    else
        error('file %s not found!',fileMetrics)
    end
%-------------------------------------------------------------------------%  
    % Standard Deviation regressor
    fileIn = fullfile(paths.EPI.dir,'2_epi.nii.gz'); %after mcflirt for this outlier detection (mask already applied)
    fileOut = fullfile(paths.EPI.epiGS,'metric_sd.txt'); 
    sentence = sprintf('%s/fslstats -t %s -S > %s',paths.FSL,fileIn,fileOut);
    [~,result] = system(sentence);
    
    SD = load(fileOut);
    
    % Zikai: when we are running multiple subjects, the configs.EPI.SDth
    % should be clear every time before we initialize scrubbing vector.
    
    % Shouldn't be fixed inside this function
    
    if isempty(configs.EPI.SDth) % if no hard threshold is defined for SD
        configs.EPI.SDth = prctile(SD,75) + (1.5*iqr(SD)); %% this follows default FSL criterion for outliers
    end
    scrubbing_sd = SD < configs.EPI.SDth;
    fprintf('-- identified %d SD-based outliers\n',nnz(~scrubbing_sd));
%-------------------------------------------------------------------------%
% Plot the scrubbed volumes
    % plot the motion metrics across volumes as bar graphs.
    figure; 
    subplot(4,1,1); plot(motionMetric_fd,'b'); 
    Subjinfo=strsplit(paths.EPI.epiGS,'/');
    title(sprintf('%s, %s, %s \n %s',Subjinfo{end-2},Subjinfo{end-1},GSreg_name,'FD'),'FontSize',8); 
    xlim([0 numTimePoints])
    subplot(4,1,2); plot(stDVARS,'r'); title('Standardized DVARS','FontSize',8); xlim([0 numTimePoints])
    subplot(4,1,3); plot(SD,'g'); title('SD','FontSize',8); xlim([0 numTimePoints])
    % plot the volumes that are scrubbed according to motion metrics.
    subplot(4,1,4)
    fdvols=find(motionMetric_fd > configs.EPI.FDth); plot(fdvols,ones(size(fdvols,1))+2,'bo'); hold on
    scvols=find(stDVARS > configs.EPI.DVARSth); plot(scvols,ones(size(scvols,1))+1,'ro'); hold on
    sdvols=find(SD > configs.EPI.SDth); plot(sdvols,ones(size(sdvols,1)),'go'); 
    title('Scrubbed Volumes by Metric','FontSize',8); xlabel('Time point (TRs)','FontSize',8);
    xlim([0 numTimePoints]); ylim([0 4])
    saveas(gcf,sprintf('%s/motion_parameters_all.png',paths.EPI.epiGS)); close all
%% -----------------------------------------------------------------------%    
    % brain mask
    volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
    % tissue masks
    volWM = MRIread(fullfile(paths.EPI.dir,'rT1_WM_mask_eroded.nii.gz'));
    volCSF = MRIread(fullfile(paths.EPI.dir,'rT1_CSF_mask_eroded.nii.gz'));
    volCSFvent = MRIread(fullfile(paths.EPI.dir,'rT1_CSFvent_mask_eroded.nii.gz'));
    % GM mask
    maskGM = MRIread(fullfile(paths.EPI.dir,'rT1_GM_mask.nii.gz'));
    volGM.vol = (maskGM.vol>0) & (volWM.vol==0) & (volCSF.vol==0) & (volBrain.vol~=0);
%-------------------------------------------------------------------------%
    % CSFvent time-series and PCA
    CSFnumVoxels = nnz(volCSFvent.vol);
    CSFmask = logical(repmat(volCSFvent.vol,[1,1,1,numTimePoints]));
    CSFts = reshape(resting.vol(CSFmask),[CSFnumVoxels,numTimePoints]); %time series of CSF voxels
    [CSFpca,CSFvar] = get_pca(CSFts',max(configs.EPI.numCompsPCA));
    CSFavg = mean(CSFts);
    CSFderiv = [0,diff(CSFavg)];
    save(fullfile(paths.EPI.epiGS,'7_dataCSF.mat'),'CSFmask','CSFts','CSFvar','CSFpca','CSFavg','CSFderiv');
    clear CSFts;
%-------------------------------------------------------------------------%
    % WM time-series and PCA
    WMnumVoxels = nnz(volWM.vol);
    WMmask = logical(repmat(volWM.vol,[1,1,1,numTimePoints]));
    WMts= reshape(resting.vol(WMmask),[WMnumVoxels,numTimePoints]); %time series of WM voxels
    [WMpca,WMvar] = get_pca(WMts',max(configs.EPI.numCompsPCA));
    WMavg = mean(WMts);
    WMderiv = [0,diff(WMavg)];
    save(fullfile(paths.EPI.epiGS,'7_dataWM.mat'),'WMmask','WMts','WMvar','WMpca','WMavg','WMderiv');
    clear WMts;
%-------------------------------------------------------------------------%
    % GM time-series and PCA
    GMnumVoxels = nnz(volGM.vol);
    GMmask = logical(repmat(volGM.vol,[1,1,1,numTimePoints]));
    GMts= reshape(resting.vol(GMmask),[GMnumVoxels,numTimePoints]); %time series of WM voxels
    [GMpca,GMvar] = get_pca(GMts',max(configs.EPI.numCompsPCA));
    GMavg = mean(GMts);
    GMderiv = [0,diff(GMavg)]; %#ok<*NASGU>
    save(fullfile(paths.EPI.epiGS,'7_dataGM.mat'),'GMmask','GMts','GMvar','GMpca','GMavg','GMderiv');
%-------------------------------------------------------------------------%
    % Global signal (GS)
    GSmask = logical(repmat(volBrain.vol,[1,1,1,numTimePoints]));
    GSnumVoxels = nnz(volBrain.vol);
    GSts = reshape(resting.vol(GSmask),[GSnumVoxels,numTimePoints]);
    [GSpca,GSvar] = get_pca(GSts',max(configs.EPI.numCompsPCA));
    GSavg = mean(GSts);
    GSderiv = [0,diff(GSavg)];
    save(fullfile(paths.EPI.epiGS,'7_dataGS.mat'),'GSmask','GSts','GSvar','GSpca','GSavg','GSderiv');
    clear GSts;

    resting.vol(~GSmask)=0;
%% -----------------------------------------------------------------------%
    % scrubbing index. 
    %   1=volume stays. 0=volume is dropped.
    scrubbing = true(numTimePoints,1);
    scrubbing = scrubbing & scrubbing_fd;
    scrubbing = scrubbing & scrubbing_dvars;
    scrubbing = scrubbing & scrubbing_sd;

    fprintf('-- identified %d motion metric-based outliers\n',nnz(~scrubbing));
    
    scrubbing(1:TR.init-1)=false;
    scrubbing((TR.end+1):end)=false;
    
    fprintf('-- identified %d start-to-end outliers \n',nnz(~scrubbing));
    
    if (configs.EPI.scrubdilate == 2)
        scrubbing = scrubbing & [scrubbing(2:end);1] & [1;scrubbing(1:end-1)]; %prev and next marked as bad too
    elseif (configs.EPI.scrubdilate == 1)
        warning('Not implemented.')
        return
    end

    fprintf('-- identified %d final outliers\n',nnz(~scrubbing));
    fprintf('-- identified %4.2f vs. %4.2f permitted outliers\n',...
        nnz(~scrubbing)/length(scrubbing),configs.EPI.scrubmaxfrac);
    
    if (nnz(scrubbing)/length(scrubbing)) >= configs.EPI.scrubmaxfrac 
        % detrend and merge regressors
        
        %regressors = [motion,motion_deriv,... % R R' regressors
        %    GSavg',GSderiv',WMavg',WMderiv',CSFavg',CSFderiv',ones(numTimePoints,1)]; % tissue related regressors
        if flags.EPI.GS==1
            regressors = [motion,motion_deriv,... % R R' regressors
            GSavg',GSderiv',WMavg',WMderiv',CSFavg',CSFderiv']; % tissue related regressors
        elseif flags.EPI.GS==0
            regressors = [motion,motion_deriv,... % R R' regressors
            WMavg',WMderiv',CSFavg',CSFderiv']; % tissue related regressors
        else
            warning('flags.EPI.GS not specified. Exiting...')
        end
        regressors_norm = nan(size(regressors));
        for i=1:size(regressors,2)
            reg = regressors(:,i);
            regMean = mean(reg(scrubbing)); %calculate mean on good points
            regStd = std(reg(scrubbing)); %calculate std on good points
            zreg = (reg-regMean)./regStd; %zscore all
            coeffs = polyfit((1:1:nnz(scrubbing))',zreg(scrubbing),1); %calculate linear trend on good points
            regressors_norm(:,i) = zreg - ((coeffs(1)*(1:length(scrubbing))) + coeffs(2))'; %detrend all
        end
        regressors(:,end+1)=1;
        regressors_norm(:,end+1)=1;
        
        save(fullfile(paths.EPI.epiGS,'7_regressors.mat'),'regressors','regressors_norm','scrubbing','scrubbing_fd','scrubbing_dvars','scrubbing_sd');
%-------------------------------------------------------------------------%
        % regress-out motion regressors  
        resting.vol = apply_regressors(resting.vol,volBrain.vol,regressors_norm,scrubbing);

        fileOut = fullfile(paths.EPI.epiGS,'7_epi.nii.gz');
        MRIwrite(resting,fileOut,'double');

        GSts_resid = reshape(resting.vol(GSmask),[GSnumVoxels,numTimePoints]); %time series of WM voxels
        GMts_resid = reshape(resting.vol(GMmask),[GMnumVoxels,numTimePoints]); %time series of WM voxels
        WMts_resid = reshape(resting.vol(WMmask),[WMnumVoxels,numTimePoints]); %time series of WM voxels
        CSFts_resid = reshape(resting.vol(CSFmask),[CSFnumVoxels,numTimePoints]); %time series of CSF voxels

        save(fullfile(paths.EPI.epiGS,'7_dataGSresid.mat'),'GSts_resid','GSmask');
        save(fullfile(paths.EPI.epiGS,'7_dataGMresid.mat'),'GMts_resid','GMmask');
        save(fullfile(paths.EPI.epiGS,'7_dataWMresid.mat'),'WMts_resid','WMmask');
        save(fullfile(paths.EPI.epiGS,'7_dataCSFresid.mat'),'CSFts_resid','CSFmask');
        save(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing','scrubbing_fd','scrubbing_dvars','scrubbing_sd');


        clear GM_ts GSts_resid GMts_resid WMts_resid CSFts_resid;
    else
        scrubmax = num2str(configs.EPI.scrubmaxfrac*100);
        sprintf('subject %s \n had >%s percent of volumes to be scrubbed! (%d out of %d volumes)',paths.subject,scrubmax,(length(scrubbing)-nnz(scrubbing)),length(scrubbing))
        if exist(fullfile(paths.EPI.epiGS,'7_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'7_epi.nii.gz'));
        end
        save(fullfile(paths.EPI.epiGS,'7_regressors.mat'),'scrubbing');
        save(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing','scrubbing_fd','scrubbing_dvars','scrubbing_sd');
    end
end

%% 8. BANDPASS
if flags.EPI.BandPass==1
    if flags.EPI.GS==1
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
    elseif flags.EPI.GS==0
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
    else
        warning('flags.EPI.GS not specified. Exiting...')
    end
    
    disp('-----------')
    disp('8. Bandpass')
    disp('-----------')

 if exist(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'file') 
    load(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing')
 else
     warning('File %s',fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'does not exist. Skipping further analysis')
        return
 end
%-------------------------------------------------------------------------% 
    if (nnz(scrubbing)/length(scrubbing)) >= configs.EPI.scrubmaxfrac    
        order = 1;
        f1 = (configs.EPI.fMin*2)*params.TR;
        f2 = (configs.EPI.fMax*2)*params.TR;
        load(fullfile(paths.EPI.epiGS,'7_dataGSresid.mat'),'GSts_resid','GSmask');
        [tsf] = apply_butterworth_filter(GSts_resid',order,f1,f2);
        clear GSts_resid;

        fileIn = fullfile(paths.EPI.epiGS,'7_epi.nii.gz');
        resting = MRIread(fileIn);
        resting.vol(GSmask) = tsf';
        fileOut = fullfile(paths.EPI.epiGS,'8_epi.nii.gz');
        err = MRIwrite(resting,fileOut,'double');
        clear GSmask tsf resting;
    else
        if exist(fullfile(paths.EPI.epiGS,'8_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'8_epi.nii.gz'));
            sprintf('Bandpass filtering for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100))
        end
    end
end

%% 9. TISSUE REGRESSORS
if flags.EPI.TissueRegressors==1
    
    disp('--------------------')
    disp('9. Tissue Regressors')
    disp('--------------------')
    
    if flags.EPI.GS==1
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
    elseif flags.EPI.GS==0
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
    else
        warning('flags.EPI.GS not specified. Exiting...')
    end
    
    
    if exist(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'file')
        load(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing')
    else
        warning('File %s',fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'does not exist. Skipping further analysis')
        return
    end
    %-------------------------------------------------------------------------%
    if (nnz(scrubbing)/length(scrubbing)) >= configs.EPI.scrubmaxfrac
        
        for nPCA = configs.EPI.numCompsPCA
            PCA_dir = sprintf('PCA%d', nPCA);
            if ~exist(fullfile(paths.EPI.epiGS, PCA_dir),'dir')
                mkdir(paths.EPI.epiGS, PCA_dir);
            end
        end
        
        resting = MRIread(fullfile(paths.EPI.epiGS,'8_epi.nii.gz'));
        zresting = resting;
        zresting.vol = zscore(zresting.vol,0,4);
        
        [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);
        %-------------------------------------------------------------------------%
        % Brain Mask
        volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
        maskBrain4D = logical(repmat(volBrain.vol,[1,1,1,numTimePoints]));
        resting.vol(~maskBrain4D)=0;
        clear maskBrain4D;
        %-------------------------------------------------------------------------%
        % tissue masks
        volGS = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
        volGM = MRIread(fullfile(paths.EPI.dir,'rT1_GM_mask.nii.gz'));
        volWM = MRIread(fullfile(paths.EPI.dir,'rT1_WM_mask_eroded.nii.gz'));
        volCSF = MRIread(fullfile(paths.EPI.dir,'rT1_CSFvent_mask_eroded.nii.gz'));
        %-------------------------------------------------------------------------%
        % CSF time-series and PCA
        CSFnumVoxels = nnz(volCSF.vol);
        CSFmask = logical(repmat(volCSF.vol,[1,1,1,numTimePoints]));
        
        CSFts = reshape(resting.vol(CSFmask),[CSFnumVoxels,numTimePoints])'; %time series of CSF voxels
        [CSFpca,CSFvar] = get_pca(CSFts,max(configs.EPI.numCompsPCA));
        CSFavg = mean(CSFts'); %#ok<*UDIM>
        save(fullfile(paths.EPI.epiGS,'9_dataCSF.mat'),'CSFmask','CSFts','CSFvar','CSFpca','CSFavg');
        %-------------------------------------------------------------------------%
        % WM time-series and PCA
        WMnumVoxels = nnz(volWM.vol);
        WMmask = logical(repmat(volWM.vol,[1,1,1,numTimePoints]));
        WMts= reshape(resting.vol(WMmask),[WMnumVoxels,numTimePoints])'; %time series of WM voxels
        [WMpca,WMvar] = get_pca(WMts,max(configs.EPI.numCompsPCA));
        WMavg = mean(WMts');
        save(fullfile(paths.EPI.epiGS,'9_dataWM.mat'),'WMmask','WMts','WMvar','WMpca','WMavg');
        %-------------------------------------------------------------------------%
        % GM time-series and PCA
        GMnumVoxels = nnz(volGM.vol);
        GMmask = logical(repmat(volGM.vol,[1,1,1,numTimePoints]));
        GMts= reshape(resting.vol(GMmask),[GMnumVoxels,numTimePoints])'; %time series of WM voxels
        [GMpca,GMvar] = get_pca(GMts,max(configs.EPI.numCompsPCA));
        GMavg = mean(GMts');
        save(fullfile(paths.EPI.epiGS,'9_dataGM.mat'),'GMmask','GMts','GMvar','GMpca','GMavg');
        %-------------------------------------------------------------------------%
        % GS time-series and PCA
        GSnumVoxels = nnz(volGS.vol);
        GSmask = logical(repmat(volGS.vol,[1,1,1,numTimePoints]));
        GSts= reshape(resting.vol(GSmask),[GSnumVoxels,numTimePoints])'; %time series of WM voxels
        [GSpca,GSvar] = get_pca(GSts,max(configs.EPI.numCompsPCA));
        GSavg = mean(GSts');
        save(fullfile(paths.EPI.epiGS,'9_dataGS.mat'),'GSmask','GSts','GSvar','GSpca','GSavg');
        %% -----------------------------------------------------------------------%
        % PCA0
        if ismember(0, configs.EPI.numCompsPCA)
            fileOut = fullfile(paths.EPI.epiGS,'PCA0','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');

            resting_vol_nopca = resting.vol;
        end
        %-------------------------------------------------------------------------%
        % PCA1 regressout tissue regressors
        if ismember(1, configs.EPI.numCompsPCA)
            if flags.EPI.GS==1
                regressors = [GSpca(:,1:1),WMpca(:,1:1),CSFpca(:,1:1),ones(numTimePoints,1)];
            elseif flags.EPI.GS==0
                regressors = [WMpca(:,1:1),CSFpca(:,1:1),ones(numTimePoints,1)];
            else
                warning('flags.EPI.GS not specified. Exiting...')
            end
            resting.vol = apply_regressors(resting_vol_nopca,volBrain.vol,regressors);
            fileOut = fullfile(paths.EPI.epiGS,'PCA1','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');
        end
        
        %-------------------------------------------------------------------------%
        % PCA2 regressout tissue regressors
        if ismember(2, configs.EPI.numCompsPCA)
            if flags.EPI.GS==1
                regressors = [GSpca(:,1:2),WMpca(:,1:2),CSFpca(:,1:2),ones(numTimePoints,1)];
            elseif flags.EPI.GS==0
                regressors = [WMpca(:,1:2),CSFpca(:,1:2),ones(numTimePoints,1)];
            else
                warning('flags.EPI.GS not specified. Exiting...')
            end
            resting.vol = apply_regressors(resting_vol_nopca,volBrain.vol,regressors);
            fileOut = fullfile(paths.EPI.epiGS,'PCA2','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');
        end
        
        %-------------------------------------------------------------------------%
        % PCA3 regressout tissue regressors
        if ismember(3, configs.EPI.numCompsPCA)
            if flags.EPI.GS==1
                regressors = [GSpca(:,1:3),WMpca(:,1:3),CSFpca(:,1:3),ones(numTimePoints,1)];
            elseif flags.EPI.GS==0
                regressors = [WMpca(:,1:3),CSFpca(:,1:3),ones(numTimePoints,1)];
            else
                warning('flags.EPI.GS not specified. Exiting...')
            end
            resting.vol = apply_regressors(resting_vol_nopca,volBrain.vol,regressors);
            fileOut = fullfile(paths.EPI.epiGS,'PCA3','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');
        end
        
        %-------------------------------------------------------------------------%
        % PCA4 regressout tissue regressors
        if ismember(4, configs.EPI.numCompsPCA)
            if flags.EPI.GS==1
                regressors = [GSpca(:,1:4),WMpca(:,1:4),CSFpca(:,1:4),ones(numTimePoints,1)];
            elseif flags.EPI.GS==0
                regressors = [WMpca(:,1:4),CSFpca(:,1:4),ones(numTimePoints,1)];
            else
                warning('flags.EPI.GS not specified. Exiting...')
            end
            resting.vol = apply_regressors(resting_vol_nopca,volBrain.vol,regressors);
            fileOut = fullfile(paths.EPI.epiGS,'PCA4','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');
        end
        
        %-------------------------------------------------------------------------%
        % PCA5 regressout tissue regressors
        if ismember(5, configs.EPI.numCompsPCA)
            if flags.EPI.GS==1
                regressors = [GSpca(:,1:5),WMpca(:,1:5),CSFpca(:,1:5),ones(numTimePoints,1)];
            elseif flags.EPI.GS==0
                regressors = [WMpca(:,1:5),CSFpca(:,1:5),ones(numTimePoints,1)];
            else
                warning('flags.EPI.GS not specified. Exiting...')
            end
            resting.vol = apply_regressors(resting_vol_nopca,volBrain.vol,regressors);
            fileOut = fullfile(paths.EPI.epiGS,'PCA5','9_epi.nii.gz');
            MRIwrite(resting,fileOut,'double');
        end
    else
        if ismember(0, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA0','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA0','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
        if ismember(1, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA1','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA1','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
        if ismember(2, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA2','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA2','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
        if ismember(3, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA3','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA3','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
        if ismember(4, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA4','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA4','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
        if ismember(5, configs.EPI.numCompsPCA) && exist(fullfile(paths.EPI.epiGS,'PCA5','9_epi.nii.gz'),'file')
            delete(fullfile(paths.EPI.epiGS,'PCA5','9_epi.nii.gz'));
            sprintf('Tissue Regressors for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
        end
    end
end

%% 9. SMOOTHING
if flags.EPI.SpatialSmooth == 1
    disp('-------------')
    disp('10. Smoothing')
    disp('-------------')
    
    if flags.EPI.GS==1
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
    elseif flags.EPI.GS==0
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
    else
        warning('flags.EPI.GS not specified. Exiting...')
    end
    
    if exist(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'file')
        load(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing')
    else
        warning('File %s',fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'does not exist. Skipping further analysis')
        return
    end
%-------------------------------------------------------------------------%  
    if (nnz(scrubbing)/length(scrubbing)) >= configs.EPI.scrubmaxfrac    
        if (configs.EPI.fwhm>0)
            sigma = configs.EPI.fwhm/(sqrt(8*log(2)));
            
            for nPCA = configs.EPI.numCompsPCA
                PCA_dir = sprintf('PCA%d', nPCA);
                fileIn = fullfile(paths.EPI.epiGS,PCA_dir,'9_epi.nii.gz');
                fileOut = fullfile(paths.EPI.epiGS,PCA_dir,'10_epi.nii.gz');
                sentence = sprintf('%s/fslmaths %s -kernel gauss %0.4f -fmean %s',paths.FSL,fileIn,sigma,fileOut);
                [~,result] = system(sentence);
            end
%-------------------------------------------------------------------------%
        else
            for nPCA = configs.EPI.numCompsPCA
                PCA_dir = sprintf('PCA%d', nPCA);
                fileIn = fullfile(paths.EPI.epiGS,PCA_dir,'8_epi.nii.gz');
                fileOut = fullfile(paths.EPI.epiGS,PCA_dir,'9_epi.nii.gz');
                copyfile(fileIn,fileOut);   
            end
        end
%-------------------------------------------------------------------------%
    else
        for nPCA = configs.EPI.numCompsPCA
            PCA_dir = sprintf('PCA%d', nPCA);
            if exist(fullfile(paths.EPI.epiGS,PCA_dir,'10_epi.nii.gz'),'file')
                delete(fullfile(paths.EPI.epiGS,PCA_dir,'10_epi.nii.gz'));
                sprintf('Smoothing for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));            
            end 
        end
    end
end

%% 10. ROIs
if flags.EPI.ROIs==1
    disp('--------')
    disp('11. ROIs')
    disp('--------')

    if flags.EPI.GS==1
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_yes');
    elseif flags.EPI.GS==0
        paths.EPI.epiGS = fullfile(paths.EPI.dir,'GSreg_no');
    else
        warning('flags.EPI.GS not specified. Exiting...')
    end
    
    if exist(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'file')
        load(fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'scrubbing')
    else
        warning('File %s',fullfile(paths.EPI.epiGS,'7_scrubbing.mat'),'does not exist.')
        return
    end
%-------------------------------------------------------------------------%
    if (nnz(scrubbing)/length(scrubbing)) >= configs.EPI.scrubmaxfrac    
        regConfigs = {};
        for i = 1:length(configs.EPI.numCompsPCA)
            regConfigs{i} = sprintf('PCA%d', configs.EPI.numCompsPCA(i));
        end
            
        volGM = MRIread(fullfile(paths.EPI.dir,'rT1_GM_mask.nii.gz')); volGM = volGM.vol;
        volWM = MRIread(fullfile(paths.EPI.dir,'rT1_WM_mask_eroded.nii.gz')); volWM = volWM.vol;
        volCSF = MRIread(fullfile(paths.EPI.dir,'rT1_CSF_mask_eroded.nii.gz')); volCSF = volCSF.vol;
        volRef = MRIread(fullfile(paths.EPI.dir,'2_epi_meanvol_mask.nii.gz'));

        for iconf = 1:length(regConfigs)
            disp(regConfigs{iconf})
            resting = MRIread(fullfile(paths.EPI.epiGS,regConfigs{iconf},'10_epi.nii.gz'));
            [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);
            
            for k=1:length(parcs.pdir)
                if parcs.pnodal(k).true==1
                    disp(parcs.plabel(k).name)
                    parcGM = MRIread(fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'_clean.nii.gz')));
                    parcGM = parcGM.vol.*(~volWM).*(~volCSF).*(volRef.vol~=0);
                    numROIs = max(parcGM(:));
                    ROIs_numVoxels = nan(numROIs,1);
                    for i_roi=1:numROIs
                        ROIs_numVoxels(i_roi,1)=nnz(parcGM==i_roi);
                    end
                    restingROIs = zeros(numROIs,numTimePoints);
                    for timePoint = 1:numTimePoints
                        aux = reshape(resting.vol(:,:,:,timePoint),[sizeX,sizeY,sizeZ]);
                        for ROI = 1:numROIs
                            voxelsROI = (parcGM==ROI);  
                            restingROIs(ROI,timePoint) = mean(aux(voxelsROI));
                        end
                        if mod(timePoint,50)==0
                            fprintf('%d out of %d\n',timePoint,numTimePoints);
                        end
                    end
                    fileOut=fullfile(paths.EPI.epiGS,regConfigs{iconf},strcat('11_epi_',parcs.plabel(k).name,'_ROIs.mat'));
                    % ROIs_numVOxels is the number of voxels belonning to each node in a partition
                    % restingROIs is the average timeseries of each region
                    save(fileOut,'restingROIs','ROIs_numVoxels');
            
                   % restingVoxels = cell(numROIs,1);
                   % for timePoint = 1:numTimePoints
                   %     aux = reshape(resting.vol(:,:,:,timePoint),[sizeX,sizeY,sizeZ]);
                   %     for ROI = 1:numROIs
                   %         voxelsROI = (parcGM==ROI);  
                   %         restingVoxels{ROI}(:,timePoint) = aux(voxelsROI);
                   %     end
                   %     if mod(timePoint,50)==0
                   %         fprintf('%d out of %d\n',timePoint,numTimePoints);
                   %     end
                   % end
                end
            end
        end
%-------------------------------------------------------------------------%        
    else
        for nPCA = configs.EPI.numCompsPCA
            PCA_dir = sprintf('PCA%d', nPCA);
            if exist(fullfile(paths.EPI.epiGS,PCA_dir,'11_epi_ROIs.mat'),'file')
                delete(fullfile(paths.EPI.epiGS,PCA_dir,'11_epi_ROIs.mat'));
                sprintf('ROIs for %s \n not done due to > %s percent of volumes being scrubbed',...
                paths.subject,num2str(configs.EPI.scrubmaxfrac*100));
            end 
        end 
    end
end
