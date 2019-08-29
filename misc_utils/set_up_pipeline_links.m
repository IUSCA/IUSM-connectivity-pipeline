% SET_UP_PIPELINE_LINKS.M

% This script is meant to help the user set up the approporiate directory 
% structure to take the data through the IUSM-connectivity-pipeline.
%
% There is an assumption that each raw subject data directory contains
% subdirectories corresponding to various MRI sequence scans.
%
% This will only create links for subject directories that do not already
% exist in data_targt directory. If you want to replace existing data, we
% can work to integrate that :) 
%
% This was initially developed and tested on SIEMENS PRISMA data collected
% at Indiana University School of Medicine.
%
% Evgeny Chumin, Indiana University, BLoomington IN, 2019
%
%%
% Provide name of structure variable (with a .name field) containing subjectList that is already read into
% matlab. You may either name it subjectList, or have it dome so in the
% line below.
%
subjectList = aging_list_prisma_crxsec;
%%
%  Provide path to data directory that contains all subject raw data
%  directories.

data_raw = '/gpfs/projects/RadImagSci/CFN/projects/mci/aging/';

%  Provide path to target directory.

data_target = '/N/dc2/projects/brainconnectomics/chumin_preproc/datadir';
%% 
%  Check status of target directory
if exist(data_target,'dir')==7
    disp('Provided target directory exists!')
    disp('Proceed with caution as to not overwrite data')
else
    disp('Creating target directory.')
    mkdir(data_target)
end
%%
diary(fullfile(data_target,'..',sprintf('link_set_up_diary_%s.log',datestr(now,'yyyymmdd'))))
disp('Raw dicom data location:')
disp(data_raw)
disp(' ')
disp('Specified target directory:')
disp(data_target)
disp(' ')
%  Loop through the provided subjectList
for i=1:length(subjectList)
  subjRAWpath=fullfile(data_raw,subjectList(i).name);  
  if exist(subjRAWpath,'dir') && ~isempty(subjRAWpath)
      sequences=dir(subjRAWpath); sequences(1:2)=[];
      subjTARGETpath=fullfile(data_target,subjectList(i).name);
      if ~exist(subjTARGETpath,'dir')
          mkdir(subjTARGETpath)
          disp('Checking sequences for:')
          disp(subjectList(i).name)
          disp(' ')
          for s=1:length(sequences)
          if sequences(s).isdir == 1
                  if strfind(sequences(s).name,'MPRAGE') > 0
                      n=1;
                  elseif strfind(sequences(s).name,'mbep2d') > 0
                      n=2;
                  elseif strfind(sequences(s).name,'DTI') > 0
                      n=3;
                  else
                      n=99;
                  end
                  
              switch n
                  case 1
                      disp('Working on T1 links to:')
                      disp(sequences(s).name)
                      dcmDir=fullfile(sequences(s).folder, sequences(s).name);
                      dcmfiles=dir(fullfile(dcmDir,'*.IMA'));
                      if length(dcmfiles) >= 176
                          T1dir=fullfile(subjTARGETpath,'T1');
                          mkdir(T1dir)
                          system(['ln -s ' dcmDir ' ' T1dir '/DICOMS'])
                          if ~isempty(dir(fullfile(T1dir,'DICOMS')))
                              disp('Link to T1 dicoms created!')
                              disp(' ')
                          else
                              warning('Link not present.')
                              disp('Need to troubleshoot, cause why the heck not!')
                              disp(' ')
                          end
                      else
                          warning('Number *.IMA dicom data found in:')
                          disp(sequences(s).name)
                          warning('is less than 176.')
                          disp('If extension is different, lets integrate it, otherwise')
                          disp('If there should be less than 176 files, lets integrate it, otherwise')
                          disp('Go find your MPRAGE data! Skipping sequence...')
                          disp(' ')
                      end
                      
                      
                  case 2
                      if strfind(sequences(s).name,'rest') > 0
                      disp('Working on EPI links to:')
                      disp(sequences(s).name)
                      dcmDir=fullfile(sequences(s).folder, sequences(s).name);
                      dcmfiles=dir(fullfile(dcmDir,'*.IMA'));
                      if length(dcmfiles) == 500
                          EPIdir=fullfile(subjTARGETpath,'EPI');
                          disp('BTW only handing single EPI sessions.')
                          disp('if you have multiples, you are going to have a problem')
                          disp('Lets work on integrating multiscan EPI link creation ;)')
                          disp(' ')
                          mkdir(EPIdir)
                          system(['ln -s ' dcmDir ' ' EPIdir '/DICOMS'])
                          if ~isempty(dir(fullfile(EPIdir,'DICOMS')))
                              disp('Link to EPI dicoms created!')
                              disp(' ')
                          else
                              warning('Link not present.')
                              disp('Need to troubleshoot, cause why the heck not!')
                              disp(' ')
                          end
                          
                          % Are there also spin echo field maps?
                          disp('Looking for *SpinEchoFieldMap* in sequence names')
                          maps=dir(fullfile(sequences(s).folder,'*SpinEchoFieldMap*'));
                          if length(maps) == 2
                              for m=1:length(maps)
                                 mappath=fullfile(maps(m).folder,maps(m).name);
                                 phase=extractBetween(maps(m).name,'Map_','_MB3');
                                 phasedir=['SEFM_' phase{1} '_DICOMS'];
                                 if m==1
                                 UNWARPdir=fullfile(subjTARGETpath,'EPI','UNWARP');
                                 mkdir(UNWARPdir)
                                 end
                                 system(['ln -s ' mappath ' ' UNWARPdir '/' phasedir])
                                 if ~isempty(dir(fullfile(UNWARPdir,phasedir)))
                                    fprintf('Link to %s created!\n', phasedir)
                                    disp(' ')
                                 else
                                    warning('Link not present.')
                                    disp('Need to troubleshoot, cause why the heck not!')
                                    disp(' ')
                                 end   
                              end                          
                          else
                              disp('Did not find it :( If you think it should be there,')
                              disp('lets work to integrage your naming convention')
                              disp(' ')
                          end
                      else
                          warning('Number *.IMA dicom data found in:')
                          disp(sequences(s).name)
                          warning('is NOT 500.')
                          disp('If extension is different, lets integrate it, otherwise')
                          disp('If there should be less (or more) than 500 files, lets integrate it, otherwise')
                          disp('Go find your EPI data! Skipping sequence...')
                          disp(' ')
                      end
                      else
                          disp(sequences(s).name)
                          disp('While it is an EPI, it is not a rest sequence.')
                          disp('If you want to add task or other EPI, lets work on that :)')
                          disp(' ')
                      end
                          
                      
                  case 3
                      disp('Working on DWI links to:')
                      disp(sequences(s).name)
                      dcmDir=fullfile(sequences(s).folder, sequences(s).name);
                      dcmfiles=dir(fullfile(dcmDir,'*.IMA'));
                      if length(dcmfiles) == 63 
                          DWIdir=fullfile(subjTARGETpath,'DWI');
                          mkdir(DWIdir)
                          system(['ln -s ' dcmDir ' ' DWIdir '/DICOMS'])
                          if ~isempty(dir(fullfile(DWIdir,'DICOMS')))
                              disp('Link to DWI dicoms created!')
                              disp(' ')
                          else
                              warning('Link not present.')
                              disp('Need to troubleshoot, cause why the heck not!')
                              disp(' ')
                          end
                          
                          % Is there a PA B0 available?
                          disp('Looking for *b0_PA in sequence names')
                          maps=dir(fullfile(sequences(s).folder,'*b0_PA'));
                          if length(maps) == 1
                              mappath=fullfile(maps(1).folder,maps(1).name);
                              phase=extractAfter(maps(1).name,'b0_');
                              phasedir=['B0_' phase '_DCM'];
                              UNWARPdir=fullfile(subjTARGETpath,'DWI','UNWARP');
                              mkdir(UNWARPdir)
                              system(['ln -s ' mappath ' ' UNWARPdir '/' phasedir])
                              if ~isempty(dir(fullfile(UNWARPdir,phasedir)))
                                  fprintf('Link to %s created!\n', phasedir)
                                  disp(' ')
                              else
                                  warning('Link not present.')
                                  disp('Need to troubleshoot, cause why the heck not!')
                                  disp(' ')
                              end
                          else                    
                              disp('Did not find it :( If you think it should be there,')
                              disp('lets work to integrage your naming convention')
                              disp(' ')
                          end
                      else
                          warning('Number *.IMA dicom data found in:')
                          disp(sequences(s).name)
                          warning('is NOT 63.')
                          disp('If extension is different, lets integrate it, otherwise')
                          disp('if there should be less (or more) than 63 files, lets integrate it, otherwise')
                          disp('go find your DWI data! Skipping sequence...')
                          disp(' ')
                      end
                      
                      
                  otherwise
                      disp(sequences(s).name)
                      disp('No match found. If this is in error,')
                      disp('let find a way to fix this :)')
                      disp(' ')
              end
          end
          end
      else
          warning('Target directory exists!')
          fprintf('Skipping %s and moving on to next subject\n',subjectList(i).name)
          disp(' ')
          if ~exist('subj_targetExists','var')
              subj_targetExists=subjectList(i);
          else
              subj_targetExists(end+1)=subjectList(i);
          end 
      end
  else
      warning('Raw data directory is nonexistent or empty!')
      fprintf('Skipping %s and moving on to next subject\n',subjectList(i).name)
      disp(' ')
      if ~exist('subj_noRAWdata','var')
          subj_noRAWdata=subjectList(i);
      else
          subj_noRAWdata(end+1)=subjectList(i);
      end
  end
end
if exist('subj_targetExists','var')
    save('link_set_up_subj_targetExists.mat','subj_targetExists')
end
if exist('subj_noRAWdata','var')
    save('link_set_up_subj_noRAWdata.mat','subj_noRAWdata')
end
diary off