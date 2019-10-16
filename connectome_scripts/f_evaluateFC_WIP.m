function [paths,flags,configs,parcs]=f_evaluateFC(paths,flags,configs,parcs,subjectinfo)
%                           F_EVALUATE_FC
% This function uses the processed EPI data to generate connectivity
% matrices as well as plot a set of figures, which should be used to
% assess quality of the image data and connectomes.
%
% Contributors:
%   Joaquin Goni, Purdue University
%   Joey Contreras, University Southern California
%   Mario Dzemidzic, Indiana University School of Medicine
%   Evgeny Chumin, Indiana University School of Medicine
%
%% Break down subjectinfo
SubjCell = strsplit(subjectinfo,',');
subjID = char(SubjCell(1));
% set colormap
set(groot, 'DefaultFigureColormap', jet)
% set resolution
dpi=200;
% output directory for figures
paths.EPI.fig = fullfile(paths.EPI.epiGS,'figures');
if ~exist(paths.EPI.fig,'dir')
    mkdir(paths.EPI.fig);
end
% Read in scrubbing data
fileIn = fullfile(paths.EPI.epiGS,'7_scrubbing.mat');
if exist(fileIn,'file')
    load(fileIn,'scrubbing');
    numTimePoints = length(scrubbing); %#ok<*NODEF>
else
    warn('File %s not found. Exiting...',fileIn)
    return
end
% Read in regressor data
if exist(fullfile(paths.EPI.epiGS,'8_epi.nii.gz'),'file')
    load(fullfile(paths.EPI.epiGS,'7_dataGM.mat'),'GMts');
    load(fullfile(paths.EPI.epiGS,'7_dataGMresid.mat'),'GMts_resid');
    load(fullfile(paths.EPI.epiGS,'7_regressors.mat')); %#ok<*LOAD>

%% FIGURE 1 Motion and Tissue Regressors
    if flags.EPI.FigsMotion==1
        figure
        subplot(3,1,1);
        plot(regressors(:,1:6)); %#ok<*IDISVAR>
        xmax = size(regressors(:,1),1);
        ylabel('mm - degrees','FontSize',10); xlim([0 xmax]);
        legend({'X','Y','Z','pitch','jaw','roll'},'FontSize',4,'Location','South','Orientation','horizontal') %#ok<*LEGINTPAR>
        legend('boxoff')
        title({sprintf('%s \n %s',subjectinfo,' Motion regressors')},'FontSize',8)

        subplot(3,1,2);
        plot(regressors(:,[1,3,5]));
        ylabel('BOLD signal','FontSize',10); xlim([0 xmax]);
        legend({'GSavg', 'WMavg','CSFavg'},'FontSize',4,'Location','South','Orientation','horizontal')
        legend('boxoff')
        title({'Tissue regressors'},'FontSize',8 )

        subplot(3,1,3);
        imagesc(regressors_norm',[-3,3]);
        ylabel('Normalized regressors','FontSize',10); 
        xlabel('Time points','FontSize',10);
        title({'Normalized regressors and their derivatives'}, 'FontSize',8)
        colormap gray

        if flags.EPI.SaveFigs==1
            figName = sprintf('epi_fig1_%s',subjID);
            save_figure(gcf,paths.EPI.fig,figName,dpi)
        end
    end

%% FIGURE 2 BOLD, its regressors, its residuals, and scrubbed volumes
    if flags.EPI.FigsMotion==1
        figure
        subplot(4,1,1);
        minY = prctile(GMts(~isnan(GMts)),2.5);
        maxY = prctile(GMts(~isnan(GMts)),97.5);

        imagesc(GMts,[minY,maxY]); title({sprintf('%s \n %s',subjectinfo,'BOLD standardized and detrended')},'FontSize',8)
        ylabel('GM voxels','FontSize',10); set(gca,'ytick',[]);
        colorbar
        colormap gray

        subplot(4,1,2)
        imagesc(regressors_norm',[-2.5,2.5]);
        ylabel('Regressors','FontSize',10); set(gca,'ytick',[]);
        title({'Normalized regressors and their derivatives'},'FontSize',8)
        colorbar
        colormap gray

        subplot(4,1,3);
        imagesc(scrubbing'); title({'Scrubbing (vertical bars)'},'FontSize',8); colorbar
        set(gca,'ytick',[]);

        subplot(4,1,4);
        minY = prctile(GMts_resid(~isnan(GMts_resid)),2.5);
        maxY = prctile(GMts_resid(~isnan(GMts_resid)),97.5);
        imagesc(GMts_resid,[minY,maxY]);
        title({'BOLD residuals after regression'},'FontSize',8);
        xlabel('Time points (TRs)','FontSize',10); 
        ylabel('GM voxels','FontSize',10); set(gca,'ytick',[]);
        colorbar
        colormap gray

        if flags.EPI.SaveFigs==1
            figName = sprintf('epi_fig2_%s',subjID);
            save_figure(gcf,paths.EPI.fig,figName,dpi)
        end
    end

%% Read in PCA data for the nodal parcellations
    for k=1:length(parcs.pdir)
        if parcs.pnodal(k).true == 1 % for every nodal parcellation
            PCA0 = load(fullfile(paths.EPI.epiGS,strcat('9_epi_',parcs.plabel(k).name,'_ROIs.mat')));
            
            validROIs = PCA0.ROIs_numVoxels>-1;   
            aux = false(length(scrubbing),1);
            aux(1:configs.EPI.Step:length(scrubbing))=true;
            scrubbing = scrubbing & aux;
            all_tp = true(length(scrubbing),1);
            
            PCA0.FC = corr(PCA0.restingROIs(validROIs,scrubbing)');
            PCA0.FCall = corr(PCA0.restingROIs(validROIs,all_tp)');
            [PCA0.FCrobust,PCA0.VProbust] = f_robust_FC(PCA0.restingROIs,validROIs,scrubbing);
            [PCA0.FCrobust_all,PCA0.VProbust_all] = f_robust_FC(PCA0.restingROIs,validROIs);
            
            

%% FIGURE 3 Mean timeseries from each ROI at each PCA component
            if (flags.EPI.FigsFC==1)
                PCAs = {};
                nPCA = length(configs.EPI.numCompsPCA);
                for i = 1:nPCA
                    variableName = sprintf('PCA%d',configs.EPI.numCompsPCA(i));
                    eval(['PCAs{' num2str(i),'} = ',variableName,';']);
                end
                
                figure
                for i = 1:nPCA
                    PCA_resting_i = PCAs{i}.restingROIs;
                    subplot(nPCA,1,i);

                    minY = floor(prctile(PCA_resting_i(~isnan(PCA_resting_i)),0.5))-5;
                    maxY = ceil(prctile(PCA_resting_i(~isnan(PCA_resting_i)),99.5))+5;
                    plot(PCA_resting_i'); axis([1 numTimePoints minY maxY]); 
                    ylabel('BOLD Signal','FontSize',10); 
                    title({sprintf('%s \n ROIs PCA%d', subjectinfo,configs.EPI.numCompsPCA(i))},'FontSize',8);
                end
                if flags.EPI.SaveFigs==1
                    figName = sprintf('epi_fig3_%s',subjID);
                    save_figure(gcf,paths.EPI.fig,figName,dpi)
                end
            end

%% FIGURE 4 Pearson correlation matrices and histograms
            if (flags.EPI.FigsFC==1)
                
                figure
                size_factor = 0.02;
                hold on;
                for i = 1:nPCA
                    h1 = subplot(2,nPCA,i);
                    p=get(h1,'pos'); p(3)=p(3)+size_factor;
                    set(h1,'pos',p);
                    imagesc(PCAs{i}.FC,[configs.EPI.minVal,configs.EPI.maxVal]);
                    axis square; xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                    if i == 1
                        title({sprintf('%s \n FC PCA%d',subjectinfo,configs.EPI.numCompsPCA(i))}, 'FontSize',8); colormap jet;
                    else
                        title({sprintf('FC PCA%d', configs.EPI.numCompsPCA(i))}, 'FontSize',8); colormap jet;
                    end
                end

                mask = triu(true(size(PCAs{1}.FC)),1);
                
                for i = 1:nPCA
                    h1 = subplot(2,nPCA,i+nPCA);
                    p=get(h1,'pos'); p(3)=p(3)+size_factor;
                    set(h1,'pos',p);
                    hist(PCAs{i}.FC(mask),configs.EPI.numBins); axis square; set(gca,'ytick',[])
                    title({sprintf('FChist PCA%d', configs.EPI.numCompsPCA(i))}, 'FontSize',8); xlabel('Pearson coeff')
                end

                if flags.EPI.SaveFigs==1
                    figName = sprintf('epi_fig4_%s',subjID);
                    save_figure(gcf,paths.EPI.fig,figName,dpi)
                end
            end

%% FIGURE 5 Grouped/Ordered connectivity matrices 
            % Look for grouping .mat files
            files=dir(fullfile(paths.MNIparcs,parcs.pdir(k).name)); files(1:2)=[];
            for e=1:length(files)
                if ~isempty(strfind(files(e).name,'.mat'))
                    % read in grouping mat file
                    Parc=load(fullfile(paths.MNIparcs,parcs.pdir(k).name,files(e).name));
                    Parc.ROIs(validROIs==0)=[];
                    [~,Parc.Order]=sort(Parc.ROIs);
                    nametmp=strsplit(files(e).name,'.');
                    Parc.title=nametmp{1}; clear nametmp
                    Parc.name=strrep(Parc.title,'_','\_');
            
                    if flags.EPI.FigsFC==1
                        figure
                        hold on;
                        
                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            imagesc(PCAs{i}.FCall(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]);
                            axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                            if i == 1
                                title({sprintf('%s \n FC PCA%d All',subjectinfo,configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                                colormap jet;
                            else
                              title({sprintf('FC PCA%d All', configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                              colormap jet; 
                            end
                        end
                        
                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i+nPCA);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            imagesc(PCAs{i}.FC(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]);
                            axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                            title({sprintf('FC PCA%d Scrubbed', configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                            colormap jet; 
                        end
                        

                        Parc.numNetworks = max(Parc.ROIs);
                        RSmask = false(size(PCAs{i}.FC));
                        for i=1:Parc.numNetworks
                            aux = (Parc.ROIs(Parc.Order)==i);
                            RSmask(aux,aux)=true;
                        end
                        
                        
                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i+2*nPCA);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            spy(RSmask); axis square; 
                            title(Parc.name,'FontSize',8);
                            xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])             
                            colormap jet; 
                        end

                        if flags.EPI.SaveFigs==1
                            figName = sprintf('epi_fig5_%s_%s',subjID,Parc.title);
                            save_figure(gcf,paths.EPI.fig,figName,dpi)
                        end
                    end
                    
%% FIGURE 6 Scrubbed Matrices
                    if flags.EPI.FigsFC==1
                        figure
                        hold on;
                        
                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            imagesc(PCAs{i}.FCrobust_all(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]);
                            axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                            if i == 1
                                title({sprintf('%s \n FC PCA%d All',subjectinfo,configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                                colormap jet;
                            else
                              title({sprintf('FC PCA%d All', configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                              colormap jet; 
                            end
                        end
                        
                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i+nPCA);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            imagesc(PCAs{i}.FCrobust(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]);
                            axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                            title({sprintf('FC PCA%d Scrubbed', configs.EPI.numCompsPCA(i))}, 'FontSize',8);
                            colormap jet; 

                        end

                        Parc.numNetworks = max(Parc.ROIs);
                        RSmask = false(size(PCAs{i}.FC));
                        for i=1:Parc.numNetworks
                            aux = (Parc.ROIs(Parc.Order)==i);
                            RSmask(aux,aux)=true;
                        end

                        for i = 1:nPCA
                            h1 = subplot(3,nPCA,i+2*nPCA);
                            p=get(h1,'pos'); p(3)=p(3)+size_factor;
                            set(h1,'pos',p);
                            spy(RSmask); axis square; 
                            title(Parc.name,'FontSize',8);
                            xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])             
                            colormap jet; 
                        end

                        if flags.EPI.SaveFigs==1
                            figName = sprintf('epi_fig6_%s_%s',subjID,Parc.title);
                            save_figure(gcf,paths.EPI.fig,figName,dpi)
                        end
                    end
                    %clear Parc
                end
            end
        end
    end

%% FIGURE 7 Cross correlation of PCA 
    if flags.EPI.FigsFC==1
        ngrid = (nPCA-1) * (nPCA-1);
        figure 
        
          
        mask(PCA0.ROIs_numVoxels<25,:) = false;
        mask(:,PCA0.ROIs_numVoxels<25) = false;
        
        comparision = fliplr(combnk(nPCA:-1:1, 2));
        for i = 1:nchoosek(nPCA,2)
            ii = comparision(i,1);
            jj = comparision(i,2);
            if i == (nPCA-1)
                plot_index = i;
            else
                plot_index = sum(0:floor(i/(nPCA-1)))+i;
            end
            subplot((nPCA-1),(nPCA-1),plot_index)
            plot(PCAs{ii}.FC(mask),PCAs{jj}.FC(mask),'.'); axis square
            r=corr(PCAs{ii}.FC(mask),PCAs{jj}.FC(mask));
            if i == 1
                title({sprintf('%s \n PCA%d vs PCA%d, r = %0.2f',subjectinfo,...
                configs.EPI.numCompsPCA(ii), configs.EPI.numCompsPCA(jj),r)}, 'FontSize',8);
            else
                title({sprintf('PCA%d vs PCA%d, r = %0.2f',configs.EPI.numCompsPCA(ii), configs.EPI.numCompsPCA(jj),...
                    r)}, 'FontSize',7);
            end
            xlabel('Pearson coeff','FontSize',6); 
            ylabel('Pearson coeff','FontSize',6);
        end

        if flags.EPI.SaveFigs==1
            figName = sprintf('epi_fig7_%s',subjID);
            save_figure(gcf,paths.EPI.fig,figName,dpi)
        end
    end
%% Save Correlation Matrices (Functional Connectomes)
    if flags.EPI.SaveMats==1
        for i = 1:nPCA 
            sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,sprintf('PCA%d', configs.EPI.numCompsPCA(i)),'FC*.mat'));
            [~,result] = system(sentence); %#ok<*ASGLU>
            FCscrubbing = PCAs{i}.FC; %#ok<*NASGU>
            FCrobust = PCAs{i}.FCrobust;
            VProbust = PCAs{i}.VProbust;
            save(fullfile(paths.EPI.epiGS,sprintf('PCA%d', configs.EPI.numCompsPCA(i)),'FCscrubbing.mat'),'FCscrubbing');
            save(fullfile(paths.EPI.epiGS,sprintf('PCA%d', configs.EPI.numCompsPCA(i)),'FCrobust.mat'),'FCrobust');
            save(fullfile(paths.EPI.epiGS,sprintf('PCA%d', configs.EPI.numCompsPCA(i)),'VProbust.mat'),'VProbust');

        end
    end
else
    if exist(fullfile(paths.EPI.epiGS,'figures'),'dir')
        PCA0=[];
        PCA1=[];
        PCA3=[];
        PCA5=[];
        sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,'figures','*.*'));
        [~,result] = system(sentence);
    end
end


%% Save Block-wise Correlation Matrices
if flags.EPI.SaveBlockMats==1
    figure
    hold on;
    for i = 1:nPCA
        h1 = subplot(2,nPCA,i);
        p=get(h1,'pos'); p(3)=p(3)+size_factor;
        set(h1,'pos',p);
        colormap(redblue)
        % Calculate block-wise mean connectivity matrices
        PCA_FC_blockWise = zeros(Parc.numNetworks, Parc.numNetworks);
        for ii = 1:Parc.numNetworks
            for jj = 1:Parc.numNetworks
                PCA_FC_blockWise(ii, jj) = nanmean(PCAs{i}.FCrobust_all(Parc.ROIs == ii, Parc.ROIs == jj),'all');
            end
        end

        imagesc(PCA_FC_blockWise,[-0.4,0.4]); axis square; xlabel('Network'); ylabel('Networks');
        set(gca,'xtick',[]); set(gca,'ytick',[])
        if i == 1
            title({sprintf('%s \n FC PCA%d All',subjectinfo, configs.EPI.numCompsPCA(i))}, 'FontSize',8);
        else
            title({sprintf('FC PCA%d All',configs.EPI.numCompsPCA(i))}, 'FontSize',8);
        end
    end

    for i = 1:nPCA
        h1 = subplot(2,nPCA,i+nPCA);
        p=get(h1,'pos'); p(3)=p(3)+size_factor;
        set(h1,'pos',p);
        % Calculate block-wise mean connectivity matrices
        PCA_FC_blockWise = zeros(Parc.numNetworks, Parc.numNetworks);
        for ii = 1:Parc.numNetworks
            for jj = 1:Parc.numNetworks
                PCA_FC_blockWise(ii, jj) = nanmean(PCAs{i}.FCrobust(Parc.ROIs == ii, Parc.ROIs == jj),'all');
            end
        end
        imagesc(PCA_FC_blockWise,[-0.4,0.4]); axis square; xlabel('Network'); ylabel('Networks');  set(gca,'xtick',[]); set(gca,'ytick',[])
        title({sprintf('FC PCA%d Scrubbed',configs.EPI.numCompsPCA(i))}, 'FontSize',8);
    end
            
    if flags.EPI.SaveFigs==1
        figName = sprintf('epi_fig8_%s_%s_blockwise',subjID,Parc.title);
        save_figure(gcf,paths.EPI.fig,figName,dpi)
    end
end
    clear Parc
end

