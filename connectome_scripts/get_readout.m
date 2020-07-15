function [RT]=get_readout(paths,dicomext)
%                               GET_READOUT
% Obtain readout time from dicom data for distortion correction in FSL EDDY.
%
%                   Effective Echo Spacing (s) = 
%        1 / (BandwidthPerPixelPhaseEncode (Hz) * MatrixSizePhase)
%
%   BandwidthPerPixelPhaseEncode -> Siemens dicom tag (0019, 1028)
%   MatrixSizePhase -> size of image in phase encode direction; usually the
% first number in the field (0051, 100b), AcquisitionMatrixText

%               Total Readout Time (FSL definition)
% Time from the center of the first echo to the center of the last
% =(actual number of phase encoding lines - 1)* effective echo spacing
%
% Actual number of phase encoding lines = 
%                       Image matrix in phase direction / GRAPPA factor

%   Evgeny Chumin, Indiana University School of Medicine, 2018
%   John West, Indiana University School of Medicine, 2018
%%
%
% JDW edit 03/16/2018 - edited function to allow passing of user defined
% DICOM file extension from system_and_sample_set_up.m
%
%%
% set paths and check for dicom directory
if isfield(paths, 'EPI')
    dicomPath=fullfile(paths.EPI.dir,'DICOMS');
else
    dicomPath=fullfile(paths.DWI.dir,'DICOMS');
end

    if length(dir(dicomPath))<= 2 %( '.' '..' are first 2 returns)
        filejson=fullfile(paths.EPI.dir,'0_epi.json');
        if exist(filejson,'file')
            jsonInfo=get_features_json(filejson,0,0);
            dim1=jsonInfo.acquisition_matrix(1);
            ees=jsonInfo.effective_echo_spacing;
            AccF=1; % need to double check this for GE
            anofel = dim1/AccF;
            RT = (anofel-1)*ees; 
        else
            warning('No files found in dicom directory. Import terminated.')
            return
        end
    else
        dicomFiles=dir(dicomPath);
        i=1; % The while loop finds the index of the first dicom file
        while contains(dicomFiles(i).name,dicomext)==0 %JDW user defined DICOM ext
            i=i+1;
        end
        % read in the dicom header
        dicomHeader=dicominfo(fullfile(dicomPath,dicomFiles(i).name));
        
        % Extract relevant information
        % matrix lines
        matrix=dicomHeader.Private_0051_100b;
        
        str_dim1 = strsplit(matrix, '*');
        dim1=str2double(str_dim1{2});
        % acceleration factor
        if isfield(dicomHeader, 'Private_0051_1011')
            AcqStr=dicomHeader.Private_0051_1011; %iPAT & MB header tag        
            Pstring=extractAfter(AcqStr,"p"); % extract everything after p
            AccF=str2double(Pstring(1)); % first return is the iPAT #
        else
            AccF=1; % No acceleration factor.
        end
        
        % bandwidth per pixel phase encode (Hz)
        bppe = dicomHeader.Private_0019_1028;
        
        % Effective Echo Spacing (s)
        ees = 1/(bppe*dim1);
        
        % Actual number of phase encoding lines
        anofel = dim1/AccF;
        
        % Total Readout Time (s)
        RT = (anofel-1)*ees;      
    end
        

