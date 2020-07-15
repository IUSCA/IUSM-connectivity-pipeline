function [paths,flags,configs,parcs,params]=f_functional_connectivity(paths,flags,configs,parcs,params)
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
    jsonFile=fullfile(paths.EPI.dir, '0_epi.json');
    [dcm_ext]=find_dcm_ext(paths.EPI.dcm);
    dicom_files=dir(fullfile(paths.EPI.dcm,sprintf('*.%s',dcm_ext)));
    if isempty(dcm_ext) && ~exist(jsonFile,'file')
        warning('No dicom (.IMA or .dcm) images found. Skipping further analysis')
        return
    elseif isempty(dcm_ext) && exist(jsonFile,'file')
        warning('No dicom images found. Using existing JSON information')
        configs.EPI.json = 1;
    end
%-------------------------------------------------------------------------%
    if configs.EPI.UseJson == 1
        disp('Using json file provided by dcm2niix to extract header information...')
        % Read Json files from dcm2niix output directory
        dcm2niix_json_loc = fullfile(paths.EPI.dir, '0_epi.json');
        dcmHeaders_json = get_features_json(dcm2niix_json_loc,  true, true);
        
        if isfield(dcmHeaders_json,'RepetitionTime')
            params.EPI.TR = dcmHeaders_json.RepetitionTime; % TR siements
        elseif isfield(dcmHeaders_json,'tr')
            params.EPI.TR = dcmHeaders_json.tr; % TR ge
        end
        fprintf(' -JSON: Repetition Time (TR): %d\n',params.EPI.TR)
        
        if isfield(dcmHeaders_json,'EchoTime')
            params.EPI.TE = dcmHeaders_json.EchoTime; % TE siemens
        elseif isfield(dcmHeaders_json,'te')
            params.EPI.TE = dcmHeaders_json.te; % TE ge
        end 
        fprintf(' -JSON: Echo Time (TE): %d\n',params.EPI.TE)
            
        if isfield(dcmHeaders_json,'FlipAngle')
            params.EPI.FlipAngle = dcmHeaders_json.FlipAngle; % Flip Angle siemens
        elseif isfield(dcmHeaders_json,'flip_angle')
            params.EPI.FlipAngle = dcmHeaders_json.flip_angle; % Flip Angle ge
        end
        fprintf(' -JSON: Flip Angle: %d\n',params.EPI.FlipAngle)
        
        if isfield(dcmHeaders_json,'EffectiveEchoSpacing')
            params.EPI.EffectiveEchoSpacing = dcmHeaders_json.EffectiveEchoSpacing; % Effective Echo Spacing siemens
        elseif isfield(dcmHeaders_json,'effective_echo_spacing')
            params.EPI.EffectiveEchoSpacing = dcmHeaders_json.effective_echo_spacing; %ge
        end
        fprintf(' -JSON: Effective Echo Spacing: %g\n',params.EPI.EffectiveEchoSpacing)
        
        if isfield(dcmHeaders_json,'BandwidthPerPixelPhaseEncode')
            params.EPI.BandwidthPerPixelPhaseEncode = dcmHeaders_json.BandwidthPerPixelPhaseEncode; % Band width Per Pixel Phase Encode
            fprintf(' -JSON: Band width Per Pixel Phase Encode: %f\n',dcmHeaders_json.BandwidthPerPixelPhaseEncode)
        else
            fprintf(' -JSON: Band width Per Pixel Phase Encode: unknown\n') %ge
        end
        
        % Extract Slicing time
        if isfield(dcmHeaders_json,'SliceTiming')
            params.EPI.slice_fractimes = dcmHeaders_json.SliceTiming; % siemens
        elseif isfield(dcmHeaders_json,'slice_timing')
            params.EPI.slice_fractimes = dcmHeaders_json.slice_timing; %ge
        end
        params.EPI.n_slice = size(params.EPI.slice_fractimes, 1); % number of slice
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
        fileInAP = fullfile(paths.EPI.SEFM,'AP.nii.gz');
        fileInPA = fullfile(paths.EPI.SEFM,'PA.nii.gz');
        if ~exist(fileInAP,'file') || ~exist(fileInPA,'file')
            fileNiiAP= 'AP';
            sentence=sprintf('rm -fr %s/%s.nii*',paths.EPI.SEFM,fileNiiAP);
            [~,result] = system(sentence); % remove any existing .nii images
            fileLog= sprintf('%s/dcm2niix_AP.log',paths.EPI.SEFM);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.scripts,fileNiiAP,paths.EPI.SEFM,paths.EPI.APdcm,fileLog);
            [~,result] = system(sentence); % import AP fieldmaps

            fileNiiPA= 'PA';
            sentence=sprintf('rm -fr %s/%s.nii*',paths.EPI.SEFM,fileNiiPA);
            [~,result] = system(sentence); % remove any existing .nii images
            fileLog= sprintf('%s/dcm2niix_PA.log',paths.EPI.SEFM);
            sentence=sprintf('%s/dcm2niix/dcm2niix -f %s -o %s -v y -x y %s > %s',paths.scripts,fileNiiPA,paths.EPI.SEFM,paths.EPI.PAdcm,fileLog);
            [~,result] = system(sentence); % import PA fieldmaps
            
            sentence=sprintf('gzip -f %s/AP.nii %s/PA.nii',paths.EPI.SEFM,paths.EPI.SEFM);
            [~,result] = system(sentence); % gzip fieldmap volumes
        end
        if exist(fileInAP,'file') && exist(fileInPA,'file')
            %Concatenate the AP then PA into single 4D image
            fileOut = fullfile(paths.EPI.SEFM,'sefield.nii.gz');
            if exist(fileOut,'file')
                sentence=sprintf('rm -fr %s',fileOut);
                [~,result] = system(sentence); % if it already exists; remove it.
            end
            sentence = sprintf('%s/fslmerge -tr %s %s %s %f',paths.FSL,fileOut,fileInAP,fileInPA,params.EPI.TR);
            [~,result]=system(sentence);
            
            % Generate an acqparams text file based on number of field maps.
            configs.EPI.SEreadOutTime = get_readout(paths,dcm_ext);
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
                % Run Topup
                    disp('Starting topup on sefield.nii.gz. This might take a while...')
                    sentence = sprintf('%s/topup --imain=%s --datain=%s --out=%s --fout=%s --iout=%s',...
                        paths.FSL,fileIn,fileParams,fileOutName,fileOutField,fileOutUnwarped);
                    [~,result]=system(sentence);
                
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
            [dcm_ext]=find_dcm_ext(paths.EPI.GREmagdcm);
            dicom_files=dir(fullfile(paths.EPI.GREmagdcm,sprintf('*.%s',dcm_ext)));
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
            warning('UNWARP DICOMS folders or nii images do not exist. Field Map correction failed.')
            return
        end
    end
end

%% 1. Slice timing correction
if flags.EPI.SliceTimingCorr==1
    
    disp('------------------------------------')
    disp('1. Slice Time Acquisition Correction')
    disp('------------------------------------')
    
  if params.EPI.TR > configs.EPI.minTR
        
    if exist(fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz');
        disp(' -Processing: 0_epi_unwarped.nii.gz')
    elseif exist(fullfile(paths.EPI.dir,'0_epi.nii.gz'),'file')
        fileIn = fullfile(paths.EPI.dir,'0_epi.nii.gz');
        disp(' -Processing: 0_epi.nii.gz')
    else
        disp('File 0_epi not found. Exiting...')
        return
    end

    sentence = sprintf('%s/fslreorient2std %s %s',paths.FSL,fileIn,fileIn);
    [~,result] = system(sentence);
    
    fileRef = fullfile(paths.EPI.dir,'slicetimes_frac.txt');
    fileOut = fullfile(paths.EPI.dir,'1_epi.nii.gz');
    load(fullfile(paths.EPI.dir,'0_param_dcm_hdr.mat'),'params')
    if params.EPI.slice_ord == 2 || configs.EPI.UseTcustom == 1 % Use custom interleave timing file
        sentence = sprintf('%s/slicetimer -i %s -o %s -r %0.4f --tcustom=%s',...
            paths.FSL,fileIn,fileOut,params.EPI.TR,fileRef);
    elseif params.EPI.slice_ord == 1 && configs.EPI.UseTcustom ~= 1 % Sequential acquisition
        if params.EPI.slice_rev == 0
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f';
        elseif params.EPI.slice_rev == 1
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --down';
        end
        sentence = sprintf(st_command,paths.FSL,fileIn,fileOut,params.EPI.TR);
    elseif params.EPI.slice_ord == 0 && configs.EPI.UseTcustom ~= 1 % Interleaved acquisition
        if params.EPI.slice_rev == 0
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --odd';
        elseif params.EPI.slice_rev == 1
            st_command = '%s/slicetimer -i %s -o %s -r %0.4f --odd --down';
        end
        sentence = sprintf(st_command,paths.FSL,fileIn,fileOut,params.EPI.TR);
    end
    [~,result] = system(sentence);
  else
      fprintf('Dataset TR is less than batch specified minTR\n')
      fprintf('Skipping slice time correction.\n')
  end
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
            fileIn = fullfile(paths.EPI.dir,'0_epi_unwarped.nii.gz');
            disp(' -Will use 0_epi_unwarped.nii.gz')
        end
    else
        fileIn = fullfile(paths.EPI.dir,'1_epi.nii.gz');
        disp(' -Will use the slice time corrected 1_epi.nii.gz as input')
    end
    
    % Compute motion outliers
    find_motion_outliers(fileIn,paths);
          
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
    fileOut = fullfile(paths.EPI.dir,'rT1_brain_dof6bbr.nii.gz');
    fileOmat = fullfile(paths.EPI.dir,'epi_2_T1_bbr.mat');
    fileWMseg = fullfile(paths.EPI.dir,'rT1_WM_mask');
    sentence = sprintf('%s/flirt -in %s -ref %s -out %s -omat %s -wmseg %s -cost bbr',...
        paths.FSL,fileIn,fileRef,fileOut,fileOmat,fileWMseg);
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

%% 5. Nuisance/Motion Parameter Regression
if flags.EPI.NuisanceReg > 0 && flags.EPI.NuisanceReg < 3
switch flags.EPI.NuisanceReg
    case 1
    disp('-----------------------')
    disp('5. ICA-AROMA: Denoising')
    disp('-----------------------')
    
    if exist(fullfile(paths.EPI.dir,'4_epi.nii.gz'),'file')
    
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
            [~,result]=system(sentence); %#ok<*ASGLU>
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
        sentence = sprintf('%s %s/ICA-AROMA/ICA_AROMA.py -in %s -out %s -mc %s -affmat %s -warp %s',...
            paths.python,paths.scripts,fileSmooth,paths.EPI.AROMAout,fileMovePar,fileMat,fileWarpField);
        [~,result]=system(sentence);
        disp(result)
        if ~exist(fullfile(paths.EPI.AROMAout,'denoised_func_data_nonaggr.nii.gz'),'file')
            warning(' -AROMA output file not found!')
            fprintf(' -Exiting... Check AROMA directory, it could be melodic ICA\n')
            fprintf('  found too many components and AROMA did not filter porperly\n')
            fprintf('  If that is the case fslfilt can be ran manually.')
        else    

            disp('   - Done.')
        end
        
        % compute percent variance removed from the data
        ICstats = dlmread(fullfile(paths.EPI.AROMAout,'melodic.ica/melodic_ICstats'));
        motionICs = dlmread(fullfile(paths.EPI.AROMAout,'classified_motion_ICs.txt'));
        for c=1:length(motionICs)
            peVar(:,c)=ICstats(motionICs(c),1);
            ptVar(:,c)=ICstats(motionICs(c),2);
        end
        peVar=sum(peVar);
        ptVar=sum(ptVar);
        fprintf('%.2f percent of explained variance in removed motion components.\n',peVar)
        fprintf('%.2f percent of total variance in removed motion components.\n',ptVar)
        clear ICstats motionICs peVar ptVar
    else
        disp('4_epi.nii.gz does not exist. Exiting...')
        return
    end
%-------------------------------------------------------------------------%    
    case 2
    disp('-----------------------------------')
    disp('5. Head Motion Parameter Regression')
    disp('-----------------------------------')
    
    if exist(fullfile(paths.EPI.dir,'4_epi.nii.gz'),'file')
    
        paths.EPI.HMP = fullfile(paths.EPI.dir,'HMPreg');
        if ~exist(paths.EPI.HMP,'dir')
            mkdir(paths.EPI.HMP)
        end
        
    % load 6 motion regressors
        load(fullfile(paths.EPI.dir,'motion.txt')); %#ok<*LOAD>
    % derivatives of 6 motion regressors
        motion_deriv = nan(size(motion)); %#ok<*NODEF>
        for i=1:size(motion,2)
            m_deriv = [0,diff(motion(:,i)')];
            motion_deriv(:,i)=m_deriv;
        end
        save(fullfile(paths.EPI.HMP,'motion12_regressors.mat'),'motion','motion_deriv')
        disp('saved motion regressors and temporal derivatives')
    % squared of motion parameters and derivatives
        if configs.EPI.numReg == 24
            motion_sq = motion.^2;
            motion_deriv_sq = motion_deriv.^2;
            save(fullfile(paths.EPI.HMP,'motion_sq_regressors.mat'),'motion_sq','motion_deriv_sq')
            disp('saved quadratics of motion and its derivatives')
        end
    else
        disp('4_epi.nii.gz does not exist. Exiting...')
        return
    end    
    otherwise
        disp('Invalid parameter selection for flags.EPI.NuisanceReg')
end

%% Physiological Regressors
if flags.EPI.PhysReg > 0 && flags.EPI.PhysReg < 3
switch flags.EPI.PhysReg
    case 1
    disp('---------')
    disp('aCompCor')
    disp('---------')
    case 2
    disp('--------------------------')
    disp('Mean CSF and WM Regression')
    disp('--------------------------')
end
    if flags.EPI.NuisanceReg==1
        fileIN=fullfile(paths.EPI.dir,'AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz');
        if exist(fileIN,'file')
            disp('Using AROMA output data')
            switch flags.EPI.PhysReg
                case 1
                    paths.EPI.PhRegDir=fullfile(paths.EPI.dir,'AROMA/aCompCorr');
                case 2
                    paths.EPI.PhRegDir=fullfile(paths.EPI.dir,'AROMA/PhysReg');
            end
        else
            disp('No AROMA output found. Cannot perform aCompCorr'); return
        end
    elseif flags.EPI.NuisanceReg==2
        fileIN=fullfile(paths.EPI.dir,'4_epi.nii.gz');
        if exist(fileIN,'file') && exist(fullfile(paths.EPI.dir,'HMPreg'),'dir')
            switch flags.EPI.PhysReg
                case 1
                    disp('Combining aCompCorr with HMP regressors')
                    paths.EPI.PhRegDir=fullfile(paths.EPI.dir,'HMPreg/aCompCorr');
                case 2
                    disp('Combining Mean CSF & WM signal with HMP regressors')
                    paths.EPI.PhRegDir=fullfile(paths.EPI.dir,'HMPreg/PhysReg');
            end
        else
            disp('Cannot find 4_epi and/or HMP directory. Cannot perform aCompCorr'); return
        end
    end
    if ~exist(paths.EPI.PhRegDir,'dir')
        mkdir(paths.EPI.PhRegDir)
    end 
    % read in data and masks
    resting = MRIread(fileIN);   
    [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);
    volCSFvent = MRIread(fullfile(paths.EPI.dir,'rT1_CSFvent_mask_eroded.nii.gz'));
    volWM = MRIread(fullfile(paths.EPI.dir,'rT1_WM_mask_eroded.nii.gz'));
    
    % read brain mask
    volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask.nii.gz'));
    volRef = MRIread(fullfile(paths.EPI.dir,'2_epi_meanvol_mask.nii.gz'));
    volBrain.vol = (volBrain.vol>0) & (volRef.vol~=0);
    fileOut = fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz');
    MRIwrite(volBrain,fileOut,'double');
    % fill holes in the brain mask, without changing FOV
    sentence = sprintf('%s/fslmaths %s -fillh %s',paths.FSL,fileOut,fileOut);
    [~,result] = system(sentence);
    
    volGS = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
%-------------------------------------------------------------------------%
    % CSFvent time-series
    CSFnumVoxels = nnz(volCSFvent.vol);
    CSFmask = logical(repmat(volCSFvent.vol,[1,1,1,numTimePoints]));
    CSFts = reshape(resting.vol(CSFmask),[CSFnumVoxels,numTimePoints]); %time series of CSF voxels
    % WM time-series
    WMnumVoxels = nnz(volWM.vol);
    WMmask = logical(repmat(volWM.vol,[1,1,1,numTimePoints]));
    WMts= reshape(resting.vol(WMmask),[WMnumVoxels,numTimePoints]); %time series of WM voxels
    % Global Signal time-series
    GSnumVoxels = nnz(volGS.vol);
    GSmask = logical(repmat(volGS.vol,[1,1,1,numTimePoints]));
    GSts= reshape(resting.vol(GSmask),[GSnumVoxels,numTimePoints]);
    
    switch flags.EPI.PhysReg
        case 1
            [CSFpca,CSFvar] = get_pca(CSFts',5);
            [WMpca,WMvar] = get_pca(WMts',5);
            save(fullfile(paths.EPI.PhRegDir,'dataPCA_WM-CSF.mat'),'CSFpca','CSFvar','CSFmask','CSFts',...
                'WMpca','WMvar','WMmask','WMts');
            disp('saved aCompCor PCA regressors')
        case 2
            CSFavg = mean(CSFts);
            CSFderiv = [0,diff(CSFavg)]';
            CSFavg = CSFavg';
            CSFavg_sq = CSFavg.^2;
            CSFderiv_sq = CSFderiv.^2;
            WMavg = mean(WMts);
            WMderiv = [0,diff(WMavg)]';
            WMavg = WMavg';
            WMavg_sq = WMavg.^2;
            WMderiv_sq = WMderiv.^2;
            save(fullfile(paths.EPI.PhRegDir,'dataMnRg_WM-CSF.mat'),'CSFavg','CSFavg_sq','CSFderiv',...
                'CSFderiv_sq','WMavg','WMavg_sq','WMderiv','WMderiv_sq');
            disp('saved mean CSF WM signal, derivatives, and quadtatics')
    end
end
%% Global Signal Regression
    if flags.EPI.GS == 1
        disp('------------------------')
        disp('Global Signal Regression')
        disp('------------------------')
        if configs.EPI.numGS > 0 && configs.EPI.numGS < 5
            GSavg = mean(GSts);
            GSderiv = [0,diff(GSavg)]';
            GSavg = GSavg';
            GSavg_sq = GSavg.^2;
            GSderiv_sq = GSderiv.^2;
            save(fullfile(paths.EPI.PhRegDir,'dataGS.mat'),'GSavg','GSavg_sq','GSderiv','GSderiv_sq')
            disp('saved global signal regressors')
        else
            disp('Invalid parameter selection for configs.EPI.numGS')
        end
    end
%% Apply Regressors
    disp('---------------------------------------------------')
    disp('Applying Regressors based on selected flags/configs')
    disp('---------------------------------------------------')
    disp('-- Creating regressor matrix with the follwing:')
    if flags.EPI.NuisanceReg == 2 
        resting = MRIread(fullfile(paths.EPI.dir,'4_epi.nii.gz'));
        volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
        switch configs.EPI.numReg
            case 24
                load(fullfile(paths.EPI.HMP,'motion12_regressors.mat'))
                load(fullfile(paths.EPI.HMP,'motion_sq_regressors.mat'))
                regressors=[motion,motion_deriv,motion_sq,motion_deriv_sq];
                disp('  -- 24 Head motion regressors')
            case 12
                load(fullfile(paths.EPI.HMP,'motion12_regressors.mat'))
                regressors=[motion,motion_deriv];
                disp('  -- 12 Head motion regressors')
        end
    elseif flags.EPI.NuisanceReg == 1
        resting = MRIread(fullfile(paths.EPI.dir,'AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz'));
        volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz'));
        regressors = double.empty;
    end
    
    if flags.EPI.GS == 1
        load(fullfile(paths.EPI.PhRegDir,'dataGS.mat'))
        switch configs.EPI.numGS
            case 4
                regressors = [regressors,GSavg,GSavg_sq,GSderiv,GSderiv_sq];
                disp('  -- 4 global signal regressors')
            case 2
                regressors = [regressors,GSavg,GSderiv];
                disp('  -- 2 global signal regressors')
            case 1
                regressors = [regressors,GSavg];
                disp('  -- 1 global signal regressor')
        end
    end
    
    if flags.EPI.PhysReg == 2
        load(fullfile(paths.EPI.PhRegDir,'dataMnRg_WM-CSF.mat'))
        switch configs.EPI.numPhys
            case 8 
                regressors = [regressors,CSFavg,CSFavg_sq,CSFderiv,...
                    CSFderiv_sq,WMavg,WMavg_sq,WMderiv,WMderiv_sq];
                disp('  -- 8 physiological regressors')
            case 4 
                regressors = [regressors,CSFavg,CSFderiv,WMavg,WMderiv];
                disp('  -- 4 physiological regressors')
            case 2 
                regressors = [regressors,CSFavg,WMavg];
                disp('  -- 2 physiological regressors')
        end
        RegressMatrix = cell.empty;
        zRegressMatrix = cell.empty;
        RegressMatrix{1,1}=regressors;
        zRegressMatrix{1,1} = zscore(RegressMatrix{1,1});
    elseif flags.EPI.PhysReg == 1
        load(fullfile(paths.EPI.PhRegDir,'dataPCA_WM-CSF.mat'))
        RegressMatrix = cell.empty;
        zRegressMatrix = cell.empty;
        disp('  -- aCompCor PC of WM & CSF regressors')
        if isempty(configs.EPI.numPC)
            disp('    -- Applying all levels of PCA removal')
            for ic = 0:5
                if ic == 0 
                    RegressMatrix{ic+1,1}=regressors;
                    zRegressMatrix{ic+1,1} = zscore(RegressMatrix{ic+1,1});
                else
                    RegressMatrix{ic+1,1}=[regressors,CSFpca(:,1:ic),WMpca(:,1:ic)];
                    zRegressMatrix{ic+1,1} = zscore(RegressMatrix{ic+1,1});
                end
            end
        elseif configs.EPI.numPC > 0 && configs.EPI.numPC < 6
            fprintf('    -- Writing prespecified removal of %d components\n',configs.EPI.numPC)
            RegressMatrix{1,1}=regressors;
            zRegressMatrix{1,1} = zscore(RegressMatrix{1,1});
        end
    end
    
    % set filename postfix for output image
    switch flags.EPI.NuisanceReg
        case 1
            nR = 'aroma';
        case 2
            nR = sprintf('hmp%d',configs.EPI.numReg);
    end
    switch flags.EPI.GS
        case 1
            nR = sprintf('%s_Gs%d',nR,configs.EPI.numGS);
    end
    switch flags.EPI.PhysReg
        case 1
            if isempty(configs.EPI.numPC)
                nR = [nR '_pca'];
            elseif configs.EPI.numPC > 0 && configs.EPI.numPC < 6
                nR = [nR '_pca' num2str(configs.EPI.numPC)];
            end
        case 2
            nR = sprintf('%s_mPhys%d',nR,configs.EPI.numPhys);
    end
    % regress-out motion/physilogical regressors  
    resting.vol(~GSmask)=0;

        for rg=1:length(zRegressMatrix)
            resid{rg,1} = apply_regressors(resting.vol,volBrain.vol,zRegressMatrix{rg,1}); 
        end
        
        resting.vol=[];
        % save data (for header info), regressors, and residuals
        save(fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output.mat',nR)),...
            'resting','volBrain','GSmask','RegressMatrix','zRegressMatrix','resid','nR','-v7.3')
        fprintf('saved regression output with %s\n',nR) 

%% 6. DEMEAN AND DETREND
if flags.EPI.DemeanDetrend == 1
    disp('---------------------')
    disp('6. Demean and Detrend')
    disp('---------------------')

    if exist(fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output.mat',nR)),'file')
        fileIn = fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output.mat',nR));
    else
        disp(' -No output found for batch defined nuisance regressed data for:')
        disp(paths.EPI.dir)
        return
    end
    
    % read data
    load(fileIn)
  for pc=1:length(resid)
    [sizeX,sizeY,sizeZ,numTimePoints] = size(resid{pc,1});
    for i=1:sizeX
        for j=1:sizeY
            for k=1:sizeZ
                if volBrain.vol(i,j,k)>0
                    TSvoxel = reshape(resid{pc,1}(i,j,k,:),[numTimePoints,1]);
                    % demean & detrend
                    TSvoxel_detrended = detrend(TSvoxel-mean(TSvoxel),'linear');
                    resid{pc,1}(i,j,k,:) = TSvoxel_detrended;
                end
            end
        end
        if (mod(i,25)==0)
            disp(i/sizeX)
        end
    end
    resid{pc,1}(~GSmask)=0;
  end
  save(fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output_dmdt.mat',nR)),...
        'resting','volBrain','GSmask','RegressMatrix','zRegressMatrix','resid','nR','-v7.3')  
end

%% 8. BANDPASS
if flags.EPI.BandPass==1
    disp('-----------')
    disp('7. Bandpass')
    disp('-----------')

    if exist(fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output_dmdt.mat',nR)),'file') 
        load(fullfile(paths.EPI.PhRegDir,sprintf('NuisanceRegression_%s_output_dmdt.mat',nR)));
        
    order = 1;
    f1 = (configs.EPI.fMin*2)*params.EPI.TR;
    f2 = (configs.EPI.fMax*2)*params.EPI.TR;
        for pc=1:length(resid)
            [~,~,~,numTimePoints] = size(resid{pc,1});
            GSts_resid = reshape(resid{pc,1}(GSmask),[nnz(volBrain.vol),numTimePoints]);
            [tsf] = apply_butterworth_filter(GSts_resid',order,f1,f2);
            clear GSts_resid;
            resting.vol=resid{pc,1};
            resting.vol(GSmask) = tsf';
            if length(resid)==1
                fileOut = fullfile(paths.EPI.PhRegDir,sprintf('7_epi_%s.nii.gz',nR));
            elseif length(resid)>1
                fileOut = fullfile(paths.EPI.PhRegDir,sprintf('7_epi_%s%d.nii.gz',nR,pc-1));
            end
            MRIwrite(resting,fileOut,'double');
            clear tsf fileOut
        end
    else
        disp('Residual timeseries from Nuisance Regressors not found. Exiting...')
        return
    end
end

%% Scrubbing
if flags.EPI.scrubbing == 1
    disp('-----------------')
    disp('Scrubbing Volumes')
    disp('-----------------')    
    
    scrub = scrub_motion_outliers(paths,configs);
    
    fileList = dir(fullfile(paths.EPI.PhRegDir,sprintf('7_epi_%s*.nii.gz',nR)));
    for fl=1:length(fileList)
        if ~isempty(strfind(fileList(fl).name,'scrubbed'))
        else
        mtype = extractBetween(fileList(fl).name,'epi_','.nii.gz');
    
        resting = MRIread(fullfile(fileList(fl).folder,fileList(fl).name));
        resting.vol=resting.vol(:,:,:,scrub);
          
        fileOut = fullfile(paths.EPI.PhRegDir,sprintf('7_epi_%s_scrubbed.nii.gz',mtype{1}));
        MRIwrite(resting,fileOut,'double');
        clear resting
        end
    end
end

%% 10. ROIs
if flags.EPI.ROIs==1
    disp('-------')
    disp('8. ROIs')
    disp('-------')
    
    volWM = MRIread(fullfile(paths.EPI.dir,'rT1_WM_mask_eroded.nii.gz')); volWM = volWM.vol;
    volCSF = MRIread(fullfile(paths.EPI.dir,'rT1_CSF_mask_eroded.nii.gz')); volCSF = volCSF.vol;
    volBrain = MRIread(fullfile(paths.EPI.dir,'rT1_brain_mask_FC.nii.gz')); volBrain = volBrain.vol;
    
    fileList = dir(fullfile(paths.EPI.PhRegDir,sprintf('7_epi_%s*.nii.gz',nR)));
    for fl=1:length(fileList)
        mtype = extractBetween(fileList(fl).name,'epi_','.nii.gz');
        paths.EPI.Mats = fullfile(paths.EPI.PhRegDir,sprintf('TimeSeries_%s',mtype{1}));
        if ~exist(paths.EPI.Mats,'dir')
            mkdir(paths.EPI.Mats)
        end
        
        resting = MRIread(fullfile(fileList(fl).folder,fileList(fl).name));
            [sizeX,sizeY,sizeZ,numTimePoints] = size(resting.vol);
            
        for k=1:length(parcs.pdir)
            if parcs.pnodal(k).true==1
                disp(parcs.plabel(k).name)
                parcGM = MRIread(fullfile(paths.EPI.dir,strcat('rT1_GM_parc_',parcs.plabel(k).name,'_clean.nii.gz')));
                parcGM = parcGM.vol.*(~volWM).*(~volCSF).*(volBrain~=0);
                numROIs = max(parcGM(:));
                ROIs_numVoxels = nan(numROIs,1);
                for i_roi=1:numROIs
                    ROIs_numVoxels(i_roi,1)=nnz(parcGM==i_roi);
                end
                restingROIs = zeros(numROIs,numTimePoints);
                ROIs_numNans = nan(numROIs,numTimePoints);
                for timePoint = 1:numTimePoints
                    aux = reshape(resting.vol(:,:,:,timePoint),[sizeX,sizeY,sizeZ]);
                    for ROI = 1:numROIs
                        voxelsROI = (parcGM==ROI);  
                        restingROIs(ROI,timePoint) = nanmean(aux(voxelsROI));
                        ROIs_numNans(ROI,timePoint)=nnz(isnan(aux(voxelsROI)));
                    end
                    if mod(timePoint,50)==0
                        fprintf('%d out of %d\n',timePoint,numTimePoints);
                    end
                end
                fileOut=fullfile(paths.EPI.Mats,strcat('8_epi_',parcs.plabel(k).name,'_ROIs.mat'));
                    % ROIs_numVOxels is the number of voxels belonging to each node in a partition
                    % restingROIs is the average timeseries of each region
                save(fileOut,'restingROIs','ROIs_numVoxels','ROIs_numNans');
            end
        end      
    end
end
end
end
