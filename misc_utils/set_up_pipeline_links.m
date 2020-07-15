% SET_UP_PIPELINE_LINKS.M

% This script is meant to help the user set up the approporiate directory 
% structure to take the data through the IUSM-connectivity-pipeline.
%
% There is an assumption that each raw subject data directory contains
% subdirectories corresponding to various MRI sequence scans.
%
% This will only create links for subject directories that do not already
% exist in data_target directory. If you want to replace existing data, we
% can work to integrate that :) 
%
% This was initially developed and tested on SIEMENS PRISMA data collected
% at Indiana University School of Medicine.
%
% Evgeny Chumin, Indiana University, Bloomington IN, 2019
%
%%
% Provide name of structure variable (with a .name field) containing subjectList that is already read into
% matlab. You may either name it subjectList, or have it dome so in the
% line below.
%
subjectList = iadc_all;
%%
%  Provide path to data directory that contains all subject raw data
%  directories.

data_raw = '/gpfs/projects/RadImagSci/CFN/projects/mci/aging/';

%  Provide path to target directory.

data_target = '/N/dc2/projects/brainconnectomics/IADC-IMAS-image-processing/datadir';
%% 
% possible sequence keywords
          mprage = {'MPRAGE'};
          rest = {'mbep2d','rest'};
            rest_tails={'_vv2bk','_scenc','_SBRef'};
          dti = {'DTI','advdiff','stejskal','diff_prod','DIFF'};
            dti_tails={'_ADC','_TRACEW','_FA','_TENSOR','_PA','_SBRef','_ColFA','_LOWB','_ColorFA','_EXP'};
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
      disp(' ')
      disp('Checking sequences for:')
      disp('--------------------------------------------------------------')
      disp(subjectList(i).name)
      disp('--------------------------------------------------------------')
      disp(' ')
      sequences = struct2cell(sequences);
      sequences = sequences(1,:)';
      for seq=1:3
        if seq == 1
          T1_id = cellfun(@(s)find(~cellfun('isempty',strfind(sequences,s))),mprage,'uni',0);
        elseif seq == 2
          REST_id = cellfun(@(s)find(~cellfun('isempty',strfind(sequences,s))),rest,'uni',0);
        elseif seq == 3
          DTI_id = cellfun(@(s)find(~cellfun('isempty',strfind(sequences,s))),dti,'uni',0);
        end
      end  
      for seq=1:3
        switch seq
          case 1
            disp('Working on Creating T1 Links')
            disp('----------------------------')
            for si = 1:length(T1_id)
              if ~isempty(T1_id{si})
                fprintf('Found %d matches for %s in %s\n\n',length(T1_id{si}),mprage{si},subjectList(i).name)    
                for si2 = 1:length(T1_id{si})     
                  s = T1_id{si}(si2);
                  disp('Working on T1 links to:')
                  disp(sequences{s})
                  dcmDir=fullfile(subjRAWpath, sequences{s});
                  [dcm_ext]=find_dcm_ext(dcmDir);
                  dcmfiles=dir(fullfile(dcmDir,sprintf('*.%s',dcm_ext)));
                  if length(dcmfiles) >= 160
                    T1dir=fullfile(subjTARGETpath,'T1');
                    if ~exist(T1dir,'dir')
                      mkdir(T1dir)
                      system(['ln -s ' dcmDir ' ' T1dir '/DICOMS'])
                      if ~isempty(dir(fullfile(T1dir,'DICOMS')))
                        disp('Link to T1 dicoms created!')
                        disp(' ')
                      else
                        warning('Link not present.')
                        disp('Need to troubleshoot, cause why the heck not!')
                        disp(' ')
                        if ~exist('subj_List_noT1link','var')
                          subj_List_noT1link = subjectList(i);
                        else
                          subj_List_noT1link(end+1) = subjectList(i);
                        end
                      end
                    else
                      fprintf(2,'T1 directory for %s exists\n',subjectList(i).name)
                      fprintf(2,'I refure to overwrite data! Fix it yourself!\n')
                      if ~exist('subj_List_multT1options','var')
                        subj_List_multT1options=subjectList(i);
                      else
                        subj_List_multT1options(end+1)=subjectList(i);
                      end
                    end
                  else
                    fprintf(2,'Number *.IMA dicom data found in:\n')
                    disp(sequences{s})
                    fprintf(2,'is less than 160.')
                    disp('If extension is different, lets integrate it, OR')
                    disp('If there should be less than 160 files, lets integrate it, otherwise')
                    disp('Go find your MPRAGE data! Skipping sequence...')
                    disp(' ')
                    if ~exist('subj_MPRAGE_issue','var')
                      subj_MPRAGE_issue = subjectList(i);
                    else
                      subj_MPRAGE_issue(end+1) = subjectList(i);
                    end  
                  end
                end  
              end
            end  
            
           case 2
            disp('Working on Creating REST Links')
            disp('------------------------------')
            REST_uID=double.empty;
            for si = 1:length(REST_id)
              REST_uID=vertcat(REST_uID,REST_id{si});
            end
            REST_uID=unique(REST_uID);
            fprintf('Found %d REST matches for %s\n\n',length(REST_uID),subjectList(i).name)
            RESTseq = sequences(REST_uID);
            r_dropSeq = cellfun(@(s)find(~cellfun('isempty',strfind(RESTseq,s))),rest_tails,'uni',0);
            restnot=double.empty;
            for ds = 1:length(r_dropSeq)
              restnot=vertcat(restnot,r_dropSeq{ds});
            end
            for si = 1:length(RESTseq)
              if sum(si==restnot)==0      
                  disp('Working on REST links to:')         
                  disp(RESTseq{si})
                  dcmDir=fullfile(subjRAWpath, RESTseq{si});
                  [dcm_ext]=find_dcm_ext(dcmDir);
                  dcmfiles=dir(fullfile(dcmDir,sprintf('*.%s',dcm_ext)));
                  if length(dcmfiles) > 160
                    EPIdir=fullfile(subjTARGETpath,'EPI');
                    disp(' ')
                    if ~exist(EPIdir,'dir')
                      mkdir(EPIdir)
                      system(['ln -s ' dcmDir ' ' EPIdir '/DICOMS'])
                      if ~isempty(dir(fullfile(EPIdir,'DICOMS')))
                        disp('Link to EPI dicoms created!')
                        disp(' ')
                      else
                        warning('Link not present.')
                        disp('Need to troubleshoot, cause why the heck not!')
                        disp(' ')
                        if ~exist('subj_List_noEPIlink','var')
                          subj_List_noEPIlink = subjectList(i);
                        else
                          subj_List_noEPIlink(end+1) = subjectList(i);
                        end
                      end
                          
                      % Are there also spin echo field maps?
                      disp('Looking for *SpinEchoFieldMap* in sequence names')
                      maps=dir(fullfile(subjRAWpath,'*SpinEchoFieldMap*'));
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
                            if ~exist('subj_List_noSEFMlink','var')
                              subj_List_noSEFMlink = subjectList(i);
                            else
                              subj_List_noSEFMlink(end+1) = subjectList(i);
                            end
                          end   
                        end                          
                      elseif length(maps) > 2
                        fprintf(2,'multiple sets of field maps found')
                        fprintf('You go sort out which one you want, okay :)')
                        if ~exist('subj_List_multSEFMdirs','var')
                          subj_List_multSEFMdirs = subjectList(i);
                        else
                          subj_List_multSEFMdirs(end+1) = subjectList(i);
                        end
                      else
                        disp('Did not find it :( If you think it should be there,')
                        disp('lets work to integrage your naming convention')
                        disp(' ')
                        if ~exist('subj_List_noSEFMdirs','var')
                          subj_List_noSEFMdirs = subjectList(i);
                        else
                          subj_List_noSEFMdirs(end+1) = subjectList(i);
                        end
                      end
                    else
                      fprintf(2,'EPI directory already exists for %s\n',subjectList(i).name)
                      fprintf('I am not in the bussiness of overwriting,\n')
                      fprintf('so go figure out why there are multiple EPI sessions.\n')
                    end
                  else
                    fprintf(2,'Number *.IMA dicom data found in:\n')
                    disp(RESTseq{si})
                    fprintf(2,'is less than 160.\n')
                    disp('If extension is different, lets integrate it, OR')
                    disp('If there should be less than 160 files, lets integrate it, otherwise')
                    disp('Go find your EPI data! Skipping sequence...')
                    disp(' ')
                    if ~exist('subj_REST_issue','var')
                      subj_REST_issue = subjectList(i);
                    else
                      subj_REST_issue(end+1) = subjectList(i);
                    end    
                  end
              end
            end
    
          case 3
            disp('Working on Creating DTI Links')
            disp('-----------------------------')
            DTI_uID=double.empty;
            for si = 1:length(DTI_id)
              DTI_uID=vertcat(DTI_uID,DTI_id{si});
            end
            DTI_uID=unique(DTI_uID);
            fprintf('Found %d DTI matches for %s\n\n',length(DTI_uID),subjectList(i).name)
            DTIseq = sequences(DTI_uID);
            dropSeq = cellfun(@(s)find(~cellfun('isempty',strfind(DTIseq,s))),dti_tails,'uni',0);
            dtinot=double.empty;
            for ds = 1:length(dropSeq)
              dtinot=vertcat(dtinot,dropSeq{ds});
            end
            for si = 1:length(DTIseq)
              if sum(si==dtinot)==0      
                  disp('Working on DWI links to:')
                  dcmDir=fullfile(subjRAWpath, DTIseq{si});
                  disp(DTIseq{si})
                  [dcm_ext]=find_dcm_ext(dcmDir);
                  dcmfiles=dir(fullfile(dcmDir,sprintf('*.%s',dcm_ext)));
                  if length(dcmfiles) >=48
                    DWIdir=fullfile(subjTARGETpath,'DWI');
                    if ~exist(DWIdir,'dir')
                      mkdir(DWIdir)
                      system(['ln -s ' dcmDir ' ' DWIdir '/DICOMS'])
                      if ~isempty(dir(fullfile(DWIdir,'DICOMS')))
                        disp('Link to DWI dicoms created!')
                        disp(' ')
                      else
                        fprintf(2,'Link not present.')
                        disp('Need to troubleshoot, cause why the heck not!')
                        disp(' ')
                        if ~exist('subj_List_noDTIlink','var')
                          subj_List_noDTIlink = subjectList(i);
                        else
                          subj_List_noDTIlink(end+1) = subjectList(i);
                        end
                      end
                          
                      % Is there a PA B0 available?
                      disp('Looking for *b0_PA in sequence names')
                      maps=dir(fullfile(subjRAWpath,'*b0_PA'));
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
                          fprintf('Link not present.')
                          disp('Need to troubleshoot, cause why the heck not!')
                          disp(' ')
                          if ~exist('subj_List_noDTIoppB0link','var')
                            subj_List_noDTIoppB0link = subjectList(i);
                          else
                            subj_List_noDTIoppB0link(end+1) = subjectList(i);
                          end
                        end
                      elseif length(maps) > 1
                        fprintf(2,'More than 1 possible opposite phase B0 exits')
                        fprintf('Go figure out which one you want and make the links manually')
                        if ~exist('subj_List_multDTIoppB0link','var')
                          subj_List_multDTIoppB0link = subjectList(i);
                        else
                          subj_List_multDTIoppB0link(end+1) = subjectList(i);
                        end  
                      else                    
                        disp('Did not find it :( If you think it should be there,')
                        disp('lets work to integrage your naming convention')
                        disp(' ')
                        if ~exist('subj_List_noDTIoppB0dir','var')
                          subj_List_DTIoppB0dir = subjectList(i);
                        else
                          subj_List_DTIoppB0dir(end+1) = subjectList(i);
                        end  
                      end
                    else
                      fprintf(2,'DWI directory already exist for %s,\n',subjectList(i).name)
                      fprintf('I refule to overwrite! Its up to you to sort out what should be done :)\n')
                    end
                  else
                    fprintf(2,'Number *.IMA dicom data found in:\n')
                    disp(DTIseq{si})
                    fprintf(2,'is < 48.\n')
                    disp(' ')
                    if ~exist('subj_DTI_issue','var')
                      subj_DTI_issue = subjectList(i);
                    else
                      subj_DTI_issue(end+1) = subjectList(i);
                    end
                  end
                
              end
            end
                                        
        end
      end
      else
        fprintf(2,'Target directory exists!\n')
        fprintf(2,'Skipping %s and moving on to next subject\n',subjectList(i).name)
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

save(fullfile(data_target,'../link_set_up_subjects_that_require_followup.mat'),'subj_*')

diary off