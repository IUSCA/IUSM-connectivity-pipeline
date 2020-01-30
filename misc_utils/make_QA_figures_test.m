%                     DATA_QUALITY_ANALYSIS_FIGURES
%                          make_QA_figures_test
%  This code generates Quality Analysis (QA) figures for data output by the
%  IUSM-connectivity-pipeline. This script will create a QA directory at
%  the level of the modality directories.
%
%  Requirements:
%       A system_and_sample_set_up.m script used to process the data, and
%       the subjectList containing IDs for which figures are to be
%       generated.
%
%  Contributors:
%           Evgeny Chumin, Indiana University, Bloomington, 2019
%%
% USER INPUT:
% system and sample set up file
clearvars
sssu = '/N/dc2/projects/brainconnectomics/IADC-IMAS-image-processing/batch_files/system_and_sample_set_up_template.m';
%
% Set to 1 the modalities you wish to create QA figures for:
section.T1brainmask = 1;
section.T1reg = 1;
section.T1parc = 0;
section.EPI = 0;
%
%%
run(sssu);
MNI = fullfile(paths.MNIparcs,'MNI_templates/MNI152_T1_1mm.nii.gz');

%%
for k=1:length(subjectList)
    if section.T1brainmask == 1
    disp(subjectList(k).name)
    paths.subject=fullfile(paths.data,subjectList(k).name); % path to subject
    paths.QAdir=fullfile(paths.subject,'QC_figures'); %output directory
    if ~exist(paths.QAdir,'dir')
        mkdir(paths.QAdir) % make output directory if it doesn't exist
    end
    
    %% 1-brain_mask_on_fov_denoised
    % Set paths and filenames
    Subj_T1=fullfile(paths.subject,configs.name.T1);
    T1fpath=fullfile(Subj_T1,'T1_fov_denoised.nii');
    maskfpath=fullfile(Subj_T1,'T1_brain_mask_filled.nii.gz');
    if isfile(T1fpath) && isfile(maskfpath)
    T1=MRIread(T1fpath);
    mask=MRIread(maskfpath);    
    % Select representative slices from T1 volume
    midslice=round(size(T1.vol,3)/2);
    slices=[midslice-30 midslice-15 midslice midslice+25 midslice+40];
    % initialize figure
    h=figure;
    h.Units='inches';
    h.Position=[1 1 20 5];
    % generate a grayscale loropmap with red as the highest intensity color
    cmap=colormap(gray(128));
    cmap(129,:)=[1 0 0];
        colormap(cmap)
    % For each representative slice
    for i=1:5
        subplot(1,5,i) % create plot in figure
        tslice=T1.vol(:,:,slices(i)); % select & display T1 slice
        h(1)=imagesc(tslice);
        hold on
        mslice=mask.vol(:,:,slices(i)); % select matching brain mask slice
        % set mask value to 1+ highest intensity in T1 slice
        mslice(mslice==1)=max(max(tslice))+1;
        h(2)=imagesc(mslice); % overlay mask
        h(2).AlphaData = 0.5; % set mask transparency
        set(gca,'Visible','off') % hide axes
        hold off
        clear tslice mslice
    end
    % Add title to figure and save as high resolution png
    sgtitle(sprintf('%s: T1 brain mask overlay',subjectList(k).name),'Interpreter','none')
    fileout = fullfile(paths.QAdir,'1-brain_mask_on_fov_denoised.png');
    count=length(dir(strcat(fileout(1:end-4),'*')));
    if count > 0
        fileout = fullfile(paths.QAdir,sprintf('1-brain_mask_on_fov_denoised_v%d.png',count+1));
    end 
    print(fileout,'-dpng','-r600')
    close all
    else
        disp('no T1_fov_denoised and/or T1_brain_mask found.')
    end
    end
    
    if section.T1reg == 1
    %% 2-T1-warped_contour_onMNI
    % read in template and subject data
    disp(subjectList(k).name)
    paths.subject=fullfile(paths.data,subjectList(k).name); % path to subject
    Subj_T1=fullfile(paths.subject,configs.name.T1);
    paths.QAdir=fullfile(paths.subject,'QC_figures'); %output directory
    if ~exist(paths.QAdir,'dir')
        mkdir(paths.QAdir) % make output directory if it doesn't exist
    end
    T1mnifile = fullfile(Subj_T1,'registration','T1_warped.nii.gz');
    if exist(T1mnifile,'file')
    T1mni=MRIread(fullfile(Subj_T1,'registration','T1_warped.nii.gz'));
    upperT1=.75*(max(max(max(T1mni.vol))));
    MNIt=MRIread(MNI);
    filename=fullfile(paths.QAdir,'2-T1_warped_contour_onMNI.gif');
    count=length(dir(strcat(filename(1:end-4),'*')));
    if count > 0
        filename = fullfile(paths.QAdir,sprintf('2-T1_warped_contour_onMNI_v%d.gif',count+1));
    end 
    % open figure
    h=figure; colormap(gray(128))
    for n=1:5:size(MNIt.vol,3) % for every 5th slice in MNI volume
        imagesc(T1mni.vol(:,:,n)); % plot MNI 
        hold all
        % overlay contour image of subject MNI space transforment T1
        contour(MNIt.vol(:,:,n),'LineWidth',1,'LineColor','r','LineStyle','-')
        set(gca,'XTickLabel',[],'YTickLabel',[])
        caxis([0 upperT1])
        title(sprintf('%s: MNI space T1 with MNI template contour overlay',subjectList(k).name),'Interpreter','none')
        drawnow
        % convert plots into iamges
        frame=getframe(h);
        im=frame2im(frame);
        [imind,cm]=rgb2ind(im,256);
        % write the gif file
        if n==1
            imwrite(imind,cm,filename,'gif','DelayTime',.2,'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',.2,'WriteMode','append')
        end
    end
    close all
    else
        fprintf('%s: No T1_warped.nii.gz found.\n',subjectList(k).name)
    end
    end
    if section.T1parc == 1
    %% 3-GM_parc_on_fov_denoised
    % get a list of parcellation files
    parcs=dir(fullfile(Subj_T1,'T1_GM_parc*'));
    % remove the dilated versions
    idx=double.empty;
    for j=1:length(parcs)
        if ~isempty(strfind(parcs(j).name,'dil'))
            idx(end+1)=j; %#ok<*SAGROW>
        end
    end
    parcs(idx)=[];
    % find max T1 intensity
    T1=MRIread(fullfile(Subj_T1,'T1_fov_denoised.nii'));
    Tmax=round(max(max(max(T1.vol))));
    % set representative slices
    midslice=round(size(T1.vol,3)/2);
    slices=[midslice-30 midslice-15 midslice midslice+25 midslice+40];
    
    % for each parcellation
    numparcs = length(parcs);
    maxIdx=double.empty;
    for p=1:numparcs
        T1p=MRIread(fullfile(Subj_T1,parcs(p).name)); % load parcellation
        maxIdx(end+1)=max(unique(T1p.vol)); % find max index value in parcellation

        % initialize figure
        if p==1
            h=figure;
            h.Units='inches';
            height = 5*numparcs; % set 5in for each parcellation
            h.Position=[1 1 20 height];
        end
        
        for n=1:5 % for each representatice slice
            numplot=n+(5*(p-1));
            subplot(numparcs,5,numplot)
            h(1)=imagesc(T1.vol(:,:,slices(n))); % plot T1 slice
            hold on
            % scale parcellation IDs to ID+twice the maximun T1 intensity
            % this ensures the color portion of the colormap is used
            mslice=T1p.vol(:,:,slices(n))+2*Tmax; 
            mslice(mslice<=2*Tmax)=0;
            h(2)=imagesc(mslice); % plot parcellation slice
            a=mslice; a(a>0)=0.7;
            h(2).AlphaData = a; % set transparency
            set(gca,'Visible','off') % hide axes
            hold off
            clear mslice
        end
    end
        % generate colormap that is a joined grayscale (low values) and
        % colors (high values); 2x size the number of nodes in parcellation.
        c2map=gray(128);
        c3map=lines(max(maxIdx));
        cpmap=vertcat(c2map,c3map);
        colormap(cpmap)
    sgtitle(sprintf('%s: GM parcs overlays',subjectList(k).name),'Interpreter','none')
    fileout = fullfile(paths.QAdir,'3-GM_parc_on_fov_denoised.png');
    count=length(dir(strcat(fileout(1:end-4),'*')));
    if count > 0
        fileout = fullfile(paths.QAdir,sprintf('3-GM_parc_on_fov_denoised_v%d.png',count+1));
    end 
    print(fileout,'-dpng','-r600')
    close all
    end
end
%%
if section.EPI ==1
for k=1:length(subjectList)
    disp(subjectList(k).name)
    paths.subject=fullfile(paths.data,subjectList(k).name); % path to subject
    paths.QAdir=fullfile(paths.subject,'QC_figures'); %output directory
    if ~exist(paths.QAdir,'dir')
        mkdir(paths.QAdir) % make output directory if it doesn't exist
    end
    
    %% 4-MCFLIRT-motion-parameters
    Subj_epi=dir(fullfile(paths.subject,[configs.name.epiFolder '*']));
    for e=1:length(Subj_epi)
        disp(Subj_epi(e).name)
        motion=dlmread(fullfile(paths.subject,Subj_epi(e).name,'motion.txt'));
        rmax = max(max(abs(motion(:,1:3))));
        h=figure('Units','inches','Position',[1 1 10 5]);
        h(1)=subplot(2,1,1);
        plot(zeros(length(motion),1),'k--')
        hold all
        plot(motion(:,1:3))
        l=rmax+(.25*rmax);
        ylim([-l l])
        title('rotation relative to mean'); legend('','x','y','z','Location','eastoutside')
        ylabel('radians')
        hold off
        
        tmax = max(max(abs(motion(:,4:6))));
        h(2)=subplot(2,1,2); %#ok<*NASGU>
        plot(zeros(length(motion),1),'k--')
        hold all
        plot(motion(:,4:6))
        l=tmax+(.25*tmax);
        ylim([-l l])
        title('translation relative to mean'); legend('','x','y','z','Location','eastoutside')
        ylabel('millimeters')
        hold off
        sgtitle(sprintf('%s: mcFLIRT motion parameters',subjectList(k).name),'Interpreter','none')
        fileout = fullfile(paths.QAdir,sprintf('4-%s-mcFLIRT_motion_parameters.png',Subj_epi(e).name));
        count=length(dir(strcat(fileout(1:end-4),'*')));
        if count > 0
            fileout = fullfile(paths.QAdir,sprintf('4-%s-mcFLIRT_motion_parameters_v%d.png',Subj_epi(e).name,count+1));
        end 
        print(fileout,'-dpng','-r600')
        close all
        %% 5-rT1_GM_mask_on_epiMeanVol
        % Set filenames/read in data
        MeanVol=MRIread(fullfile(paths.subject,Subj_epi(e).name,'2_epi_meanvol.nii.gz'));
        mask=MRIread(fullfile(paths.subject,Subj_epi(e).name,'rT1_GM_mask.nii.gz'));
        % Select representative slices from EPI volume
        midslice=round(size(MeanVol.vol,3)/2);
        slices=[midslice-15 midslice-5 midslice midslice+5 midslice+15];
        % initialize figure
        h=figure;
        h.Units='inches';
        h.Position=[1 1 20 5];
        % generate a grayscale loropmap with red as the highest intensity color
        cmap=colormap(gray(128));
        cmap(129,:)=[1 0 0];
            colormap(cmap)
        % For each representative slice
        for i=1:5
            subplot(1,5,i) % create plot in figure
            vslice=MeanVol.vol(:,:,slices(i)); % select & display epi slice
            h(1)=imagesc(vslice);
            hold on
            mslice=mask.vol(:,:,slices(i)); % select matching brain mask slice
            % set mask value to 1+ highest intensity in epi slice
            mslice(mslice==1)=max(max(vslice))+1;
            h(2)=imagesc(mslice); % overlay mask
            h(2).AlphaData = 0.5; % set mask transparency
            set(gca,'Visible','off') % hide axes
            hold off
            clear vslice mslice
        end
        % Add title to figure and save as high resolution png
        sgtitle(sprintf('%s: rT1_GM_mask on epi_meanvol',subjectList(k).name),'Interpreter','none')
        fileout = fullfile(paths.QAdir,sprintf('5-%s-rT1_GM_mask_on_epiMeanVol.png',Subj_epi(e).name));
        count=length(dir(strcat(fileout(1:end-4),'*')));
        if count > 0
            fileout = fullfile(paths.QAdir,sprintf('5-%s-rT1_GM_mask_on_epiMeanVol_v%d.png',Subj_epi(e).name,count+1));
        end 
        print(fileout,'-dpng','-r600')
        close all
        
        %% 6-rT1_GM_parc_on_epi_meanVol
        % get a list of parcellation files
        parcs=dir(fullfile(paths.subject,Subj_epi(e).name,'rT1_GM_parc*clean*'));
        % find max T1 intensity
        Emax=round(max(max(max(MeanVol.vol))));
        % for each parcellation
        numparcs = length(parcs);
        maxIdx=double.empty;
        for p=1:numparcs
            EPIp=MRIread(fullfile(paths.subject,Subj_epi(e).name,parcs(p).name)); % load parcellation
            maxIdx(end+1)=max(unique(EPIp.vol)); % find max index value in parcellation
            
            % initialize figure
            if p==1
                h=figure;
                h.Units='inches';
                height = 5*numparcs; % set 5in for each parcellation
                h.Position=[1 1 20 height];
            end
        
            for n=1:5 % for each representatice slice
                numplot=n+(5*(p-1));
                subplot(numparcs,5,numplot)
                h(1)=imagesc(MeanVol.vol(:,:,slices(n))); % plot T1 slice
                hold on
                % scale parcellation IDs to ID+twice the maximun T1 intensity
                % this ensures the color portion of the colormap is used
                mslice=EPIp.vol(:,:,slices(n))+Emax; 
                mslice(mslice<=Emax)=0;
                h(2)=imagesc(mslice); % plot parcellation slice
                a=mslice; a(a>0)=0.7;
                h(2).AlphaData = a; % set transparency
                set(gca,'Visible','off') % hide axes
                hold off
                clear mslice
            end
        end
        % generate colormap that is a joined grayscale (low values) and
        % colors (high values); 2x size the number of nodes in parcellation.
        c2map=gray(10000);
        c3map=lines(max(maxIdx));
        cpmap=vertcat(c2map,c3map);
        colormap(cpmap)
        sgtitle(sprintf('%s:%s: GM parc overlays',subjectList(k).name,Subj_epi(e).name),'Interpreter','none')
        fileout = fullfile(paths.QAdir,sprintf('6-%s-rT1_GM_parc_on_epi_meanVol.png',Subj_epi(e).name));
        count=length(dir(strcat(fileout(1:end-4),'*')));
        if count > 0
            fileout = fullfile(paths.QAdir,sprintf('6-%s-rT1_GM_parc_on_epi_meanVol_v%d.png',Subj_epi(e).name,count+1));
        end 
        print(fileout,'-dpng','-r600')
        close all 
        
        %% 7-ICA_AROMA_component_assessment
        figAROMA = fullfile(paths.subject,Subj_epi(e).name,'AROMA/AROMA-output/ICA_AROMA_component_assessment.pdf');
        if exist(figAROMA,'file')
            fileout = fullfile(paths.QAdir,sprintf('7-%s-ICA_AROMA_component_assessment.pdf',Subj_epi(e).name));
            count=length(dir(strcat(fileout(1:end-4),'*')));
            if count > 0
                fileout = fullfile(paths.QAdir,sprintf('7-%s-ICA_AROMA_component_assessment_v%d.pdf',Subj_epi(e).name,count+1));
            end 
            copyfile(figAROMA,fileout)
            close all 
        else
            textout=fullfile(paths.QAdir,sprintf('7-%s-no_aroma_figure_found',Subj_epi(e).name));
            system('touch %s',textout)
        end
        
        %% 8-Nuisance_regressors
        %TBD
    end
end 
end
%%

% DWI
clear subjectList












