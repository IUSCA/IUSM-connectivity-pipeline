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
fileIn = fullfile(paths.EPI.epiGS,'6_scrubbing.mat');
if exist(fileIn,'file')
    load(fileIn,'scrubbing');
    numTimePoints = length(scrubbing); %#ok<*NODEF>
else
    warn('File %s not found. Exiting...',fileIn)
    return
end
% Read in regressor data
if exist(fullfile(paths.EPI.epiGS,'PCA5','9_epi.nii.gz'),'file')
    load(fullfile(paths.EPI.epiGS,'6_dataGM.mat'),'GMts');
    load(fullfile(paths.EPI.epiGS,'6_dataGMresid.mat'),'GMts_resid');
    load(fullfile(paths.EPI.epiGS,'6_regressors.mat')); %#ok<*LOAD>

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
        plot(regressors(:,[13,15,17]));
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
            PCA0 = load(fullfile(paths.EPI.epiGS,'PCA0',strcat('10_epi_',parcs.plabel(k).name,'_ROIs.mat')));
            
            validROIs = PCA0.ROIs_numVoxels>-1;   
            aux = false(length(scrubbing),1);
            aux(1:configs.EPI.Step:length(scrubbing))=true;
            scrubbing = scrubbing & aux;
    
            PCA0.FC = corr(PCA0.restingROIs(validROIs,scrubbing)');
            PCA0.FCall = corr(PCA0.restingROIs(validROIs,:)');
            [PCA0.FCrobust,PCA0.VProbust] = f_robust_FC(PCA0.restingROIs,validROIs,scrubbing);
            [PCA0.FCrobust_all,PCA0.VProbust_all] = f_robust_FC(PCA0.restingROIs,validROIs);

            PCA1 = load(fullfile(paths.EPI.epiGS,'PCA1',strcat('10_epi_',parcs.plabel(k).name,'_ROIs.mat')));
            PCA1.FC = corr(PCA1.restingROIs(validROIs,scrubbing)');
            PCA1.FCall = corr(PCA1.restingROIs(validROIs,:)');
            [PCA1.FCrobust,PCA1.VProbust] = f_robust_FC(PCA1.restingROIs,validROIs,scrubbing);
            [PCA1.FCrobust_all,PCA1.VProbust_all] = f_robust_FC(PCA1.restingROIs,validROIs);

            PCA3 = load(fullfile(paths.EPI.epiGS,'PCA3',strcat('10_epi_',parcs.plabel(k).name,'_ROIs.mat')));
            PCA3.FC = corr(PCA3.restingROIs(validROIs,scrubbing)');
            PCA3.FCall = corr(PCA3.restingROIs(validROIs,:)');
            [PCA3.FCrobust,PCA3.VProbust] = f_robust_FC(PCA3.restingROIs,validROIs,scrubbing);
            [PCA3.FCrobust_all,PCA3.VProbust_all] = f_robust_FC(PCA3.restingROIs,validROIs);

            PCA5 = load(fullfile(paths.EPI.epiGS,'PCA5',strcat('10_epi_',parcs.plabel(k).name,'_ROIs.mat')));
            PCA5.FC = corr(PCA5.restingROIs(validROIs,scrubbing)');
            PCA5.FCall = corr(PCA5.restingROIs(validROIs,:)');
            [PCA5.FCrobust,PCA5.VProbust] = f_robust_FC(PCA5.restingROIs,validROIs,scrubbing);
            [PCA5.FCrobust_all,PCA5.VProbust_all] = f_robust_FC(PCA5.restingROIs,validROIs);

%% FIGURE 3 Mean timeseries from each ROI at each PCA component
            if (flags.EPI.FigsFC==1)
                figure
                subplot(4,1,1);
                minY = floor(prctile(PCA0.restingROIs(~isnan(PCA0.restingROIs)),0.5))-5;
                maxY = ceil(prctile(PCA0.restingROIs(~isnan(PCA0.restingROIs)),99.5))+5;
                plot(PCA0.restingROIs'); axis([1 numTimePoints minY maxY]); 
                ylabel('BOLD Signal','FontSize',10); 
                title({sprintf('%s \n %s', subjectinfo,'ROIs PCA0')},'FontSize',8);
                subplot(4,1,2);
                minY = floor(prctile(PCA1.restingROIs(~isnan(PCA1.restingROIs)),0.5))-5;
                maxY = ceil(prctile(PCA1.restingROIs(~isnan(PCA1.restingROIs)),99.5))+5;
                plot(PCA1.restingROIs'); axis([1 numTimePoints minY maxY]); 
                ylabel('BOLD Signal','FontSize',10); title({'ROIs PCA1'}, 'FontSize',8);
                subplot(4,1,3);
                minY = floor(prctile(PCA3.restingROIs(~isnan(PCA3.restingROIs)),0.5))-5;
                maxY = ceil(prctile(PCA3.restingROIs(~isnan(PCA3.restingROIs)),99.5))+5;
                plot(PCA3.restingROIs'); axis([1 numTimePoints minY maxY]); 
                ylabel('BOLD Signal','FontSize',10); title({'ROIs PCA3'}, 'FontSize',8);
                subplot(4,1,4);
                minY = floor(prctile(PCA5.restingROIs(~isnan(PCA5.restingROIs)),0.5))-5;
                maxY = ceil(prctile(PCA5.restingROIs(~isnan(PCA5.restingROIs)),99.5))+5;
                plot(PCA5.restingROIs'); axis([1 numTimePoints minY maxY]); 
                xlabel('Time points (TRs)','FontSize',10); 
                ylabel('BOLD Signal','FontSize',10); title({'ROIs PCA5'}, 'FontSize',8);

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
                h1 = subplot(2,4,1);
                p=get(h1,'pos'); p(3)=p(3)+size_factor;
                set(h1,'pos',p);
                imagesc(PCA0.FC,[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])

                title({sprintf('%s \n %s',subjectinfo,'FC PCA0')}, 'FontSize',8); colormap jet;
                h2 = subplot(2,4,2);
                p=get(h2,'pos'); p(3)=p(3)+size_factor;
                set(h2,'pos',p);
                imagesc(PCA1.FC,[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                title({'FC PCA1'}, 'FontSize',8); colormap jet;
                h3 = subplot(2,4,3);
                p=get(h3,'pos'); p(3)=p(3)+size_factor;
                set(h3,'pos',p);
                imagesc(PCA3.FC,[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                title({'FC PCA3'}, 'FontSize',8); colormap jet;
                h4 = subplot(2,4,4);
                p=get(h4,'pos'); p(3)=p(3)+size_factor;
                set(h4,'pos',p);
                imagesc(PCA5.FC,[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                title({'FC PCA5'}, 'FontSize',8); colormap jet;

                mask = triu(true(size(PCA0.FC)),1);

                h5 = subplot(2,4,5);
                p=get(h5,'pos'); p(3)=p(3)+size_factor;
                set(h5,'pos',p);
                hist(PCA0.FC(mask),configs.EPI.numBins); axis square; set(gca,'ytick',[])
                title({'FChist PCA0'}, 'FontSize',8); xlabel('Pearson coeff')
                h6 = subplot(2,4,6);
                p=get(h6,'pos'); p(3)=p(3)+size_factor;
                set(h6,'pos',p);
                hist(PCA1.FC(mask),configs.EPI.numBins); axis square; set(gca,'ytick',[])
                title({'FChist PCA1'}, 'FontSize',8); xlabel('Pearson coeff')
                h7 = subplot(2,4,7);
                p=get(h7,'pos'); p(3)=p(3)+size_factor;
                set(h7,'pos',p);
                hist(PCA3.FC(mask),configs.EPI.numBins); axis square; set(gca,'ytick',[])
                title({'FChist PCA3'}, 'FontSize',8); xlabel('Pearson coeff')
                h8 = subplot(2,4,8);
                p=get(h8,'pos'); p(3)=p(3)+size_factor;
                set(h8,'pos',p);
                hist(PCA5.FC(mask),configs.EPI.numBins); axis square; set(gca,'ytick',[])
                title({'FChist PCA5'}, 'FontSize',8); xlabel('Pearson coeff')

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
                        h1 = subplot(3,4,1);
                        p=get(h1,'pos'); p(3)=p(3)+size_factor;
                        set(h1,'pos',p);
                        imagesc(PCA0.FCall(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({sprintf('%s \n %s',subjectinfo,'FC PCA0 All')}, 'FontSize',8); colormap jet;
                        h2 = subplot(3,4,2);
                        p=get(h2,'pos'); p(3)=p(3)+size_factor;
                        set(h2,'pos',p);
                        imagesc(PCA1.FCall(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA1 All'}, 'FontSize',8); colormap jet;
                        h3 = subplot(3,4,3);
                        p=get(h3,'pos'); p(3)=p(3)+size_factor;
                        set(h3,'pos',p);
                        imagesc(PCA3.FCall(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA3 All'}, 'FontSize',8); colormap jet;
                        h4 = subplot(3,4,4);
                        p=get(h4,'pos'); p(3)=p(3)+size_factor;
                        set(h4,'pos',p);
                        imagesc(PCA5.FCall(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA5 All'}, 'FontSize',8); colormap jet;

                        h1 = subplot(3,4,5);
                        p=get(h1,'pos'); p(3)=p(3)+size_factor;
                        set(h1,'pos',p);
                        imagesc(PCA0.FC(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA0 Scrubbed'}, 'FontSize',8); colormap jet;
                        h2 = subplot(3,4,6);
                        p=get(h2,'pos'); p(3)=p(3)+size_factor;
                        set(h2,'pos',p);
                        imagesc(PCA1.FC(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA1 Scrubbed'},'FontSize',8); colormap jet;
                        h3 = subplot(3,4,7);
                        p=get(h3,'pos'); p(3)=p(3)+size_factor;
                        set(h3,'pos',p);
                        imagesc(PCA3.FC(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA3 Scrubbed'}, 'FontSize',8); colormap jet;
                        h4 = subplot(3,4,8);
                        p=get(h4,'pos'); p(3)=p(3)+size_factor;
                        set(h4,'pos',p);
                        imagesc(PCA5.FC(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA5 Scrubbed'}, 'FontSize',8); colormap jet;

                        Parc.numNetworks = max(Parc.ROIs);
                        RSmask = false(size(PCA0.FC));
                        for i=1:Parc.numNetworks
                            aux = (Parc.ROIs(Parc.Order)==i);
                            RSmask(aux,aux)=true;
                        end

                        h5 = subplot(3,4,9);
                        p=get(h5,'pos'); p(3)=p(3)+size_factor;
                        set(h5,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name,'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h6 = subplot(3,4,10);
                        p=get(h6,'pos'); p(3)=p(3)+size_factor;
                        set(h6,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name,'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h7 = subplot(3,4,11);
                        p=get(h7,'pos'); p(3)=p(3)+size_factor;
                        set(h7,'pos',p);
                        spy(RSmask); axis square;  
                        title(Parc.name, 'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h8 = subplot(3,4,12);
                        p=get(h8,'pos'); p(3)=p(3)+size_factor;
                        set(h8,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name, 'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])

                        if flags.EPI.SaveFigs==1
                            figName = sprintf('epi_fig5_%s_%s',subjID,Parc.title);
                            save_figure(gcf,paths.EPI.fig,figName,dpi)
                        end
                    end
%% FIGURE 6 Scrubbed Matrices
                    if flags.EPI.FigsFC==1
                        figure
                        hold on;
                        h1 = subplot(3,4,1);
                        p=get(h1,'pos'); p(3)=p(3)+size_factor;
                        set(h1,'pos',p);
                        imagesc(PCA0.FCrobust_all(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({sprintf('%s \n %s',subjectinfo,'FC PCA0 All')}, 'FontSize',8);
                        h2 = subplot(3,4,2);
                        p=get(h2,'pos'); p(3)=p(3)+size_factor;
                        set(h2,'pos',p);
                        imagesc(PCA1.FCrobust_all(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA1 All'}, 'FontSize',8);
                        h3 = subplot(3,4,3);
                        p=get(h3,'pos'); p(3)=p(3)+size_factor;
                        set(h3,'pos',p);
                        imagesc(PCA3.FCrobust_all(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA3 All'}, 'FontSize',8);
                        h4 = subplot(3,4,4);
                        p=get(h4,'pos'); p(3)=p(3)+size_factor;
                        set(h4,'pos',p);
                        imagesc(PCA5.FCrobust_all(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA5 All'}, 'FontSize',8);

                        h1 = subplot(3,4,5);
                        p=get(h1,'pos'); p(3)=p(3)+size_factor;
                        set(h1,'pos',p);
                        imagesc(PCA0.FCrobust(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA0 Scrubbed'}, 'FontSize',8);
                        h2 = subplot(3,4,6);
                        p=get(h2,'pos'); p(3)=p(3)+size_factor;
                        set(h2,'pos',p);
                        imagesc(PCA1.FCrobust(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA1 Scrubbed'}, 'FontSize',8);
                        h3 = subplot(3,4,7);
                        p=get(h3,'pos'); p(3)=p(3)+size_factor;
                        set(h3,'pos',p);
                        imagesc(PCA3.FCrobust(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA3 Scrubbed'}, 'FontSize',8);
                        h4 = subplot(3,4,8);
                        p=get(h4,'pos'); p(3)=p(3)+size_factor;
                        set(h4,'pos',p);
                        imagesc(PCA5.FCrobust(Parc.Order,Parc.Order),[configs.EPI.minVal,configs.EPI.maxVal]); axis square; xlabel('ROIs'); ylabel('ROIs');  set(gca,'xtick',[]); set(gca,'ytick',[])
                        title({'FC PCA5 Scrubbed'}, 'FontSize',8);

                        Parc.numNetworks = max(Parc.ROIs);
                        RSmask = false(size(PCA0.FC));
                        for i=1:Parc.numNetworks
                            aux = (Parc.ROIs(Parc.Order)==i);
                            RSmask(aux,aux)=true;
                        end

                        h5 = subplot(3,4,9);
                        p=get(h5,'pos'); p(3)=p(3)+size_factor;
                        set(h5,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name,'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h6 = subplot(3,4,10);
                        p=get(h6,'pos'); p(3)=p(3)+size_factor;
                        set(h6,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name, 'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h7 = subplot(3,4,11);
                        p=get(h7,'pos'); p(3)=p(3)+size_factor;
                        set(h7,'pos',p);
                        spy(RSmask); axis square;  
                        title(Parc.name,'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])
                        h8 = subplot(3,4,12);
                        p=get(h8,'pos'); p(3)=p(3)+size_factor;
                        set(h8,'pos',p);
                        spy(RSmask); axis square; 
                        title(Parc.name,'FontSize',8); xlabel('ROIs'); ylabel('ROIs'); set(gca,'xtick',[]); set(gca,'ytick',[])

                        if flags.EPI.SaveFigs==1
                            figName = sprintf('epi_fig6_%s_%s',subjID,Parc.title);
                            save_figure(gcf,paths.EPI.fig,figName,dpi)
                        end
                    end
                    clear Parc
                end
            end
        end
    end

%% FIGURE 7 Cross correlation of PCA 
    if flags.EPI.FigsFC==1
        figure 
        mask(PCA0.ROIs_numVoxels<25,:) = false;
        mask(:,PCA0.ROIs_numVoxels<25) = false;
        subplot(3,3,1)
        plot(PCA0.FC(mask),PCA1.FC(mask),'.'); axis square
        r=corr(PCA0.FC(mask),PCA1.FC(mask));
        title({sprintf('%s \n PCA0 vs PCA1, r = %0.2f',subjectinfo, r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);
        subplot(3,3,2);
        plot(PCA0.FC(mask),PCA3.FC(mask),'.'); axis square;
        r=corr(PCA0.FC(mask),PCA3.FC(mask));
        title({sprintf('PCA0 vs PCA3, r = %0.2f',r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);
        subplot(3,3,3);
        plot(PCA0.FC(mask),PCA5.FC(mask),'.'); axis square;
        r=corr(PCA0.FC(mask),PCA5.FC(mask));
        title({sprintf('PCA0 vs PCA5, r = %0.2f',r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);

        subplot(3,3,5);
        plot(PCA1.FC(mask),PCA3.FC(mask),'.'); axis square;
        r=corr(PCA1.FC(mask),PCA3.FC(mask));
        title({sprintf('PCA1 vs PCA3, r = %0.2f',r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);
        subplot(3,3,6);
        plot(PCA1.FC(mask),PCA5.FC(mask),'.'); axis square;
        r=corr(PCA1.FC(mask),PCA5.FC(mask));
        title({sprintf('PCA1 vs PCA5, r = %0.2f',r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);

        subplot(3,3,9);
        plot(PCA3.FC(mask),PCA5.FC(mask),'.'); axis square;
        r=corr(PCA3.FC(mask),PCA5.FC(mask));
        title({sprintf('PCA3 vs PCA5, r = %0.2f',r)}, 'FontSize',8);
        xlabel('Pearson coeff','FontSize',8); 
        ylabel('Pearson coeff','FontSize',8);

        if flags.EPI.SaveFigs==1
            figName = sprintf('epi_fig7_%s',subjID);
            save_figure(gcf,paths.EPI.fig,figName,dpi)
        end
    end
%% Save Correlation Matrices (Functional Connectomes)
    if flags.EPI.SaveMats==1
        
        sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,'PCA0','FC*.mat'));
        [~,result] = system(sentence); %#ok<*ASGLU>
        FCscrubbing = PCA0.FC; %#ok<*NASGU>
        FCrobust = PCA0.FCrobust;
        VProbust = PCA0.VProbust;
        save(fullfile(paths.EPI.epiGS,'PCA0','FCscrubbing.mat'),'FCscrubbing');
        save(fullfile(paths.EPI.epiGS,'PCA0','FCrobust.mat'),'FCrobust');
        save(fullfile(paths.EPI.epiGS,'PCA0','VProbust.mat'),'VProbust');

        sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,'PCA1','FC*.mat'));
        [~,result] = system(sentence);
        FCscrubbing = PCA1.FC;
        FCrobust = PCA1.FCrobust;
        VProbust = PCA1.VProbust;
        save(fullfile(paths.EPI.epiGS,'PCA1','FCscrubbing.mat'),'FCscrubbing');
        save(fullfile(paths.EPI.epiGS,'PCA1','FCrobust.mat'),'FCrobust');
        save(fullfile(paths.EPI.epiGS,'PCA1','VProbust.mat'),'VProbust');

        sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,'PCA3','FC*.mat'));
        [status,result] = system(sentence);
        FCscrubbing = PCA3.FC;
        FCrobust = PCA3.FCrobust;
        VProbust = PCA3.VProbust;
        save(fullfile(paths.EPI.epiGS,'PCA3','FCscrubbing.mat'),'FCscrubbing');
        save(fullfile(paths.EPI.epiGS,'PCA3','FCrobust.mat'),'FCrobust');
        save(fullfile(paths.EPI.epiGS,'PCA3','VProbust.mat'),'VProbust');

        sentence = sprintf('rm %s',fullfile(paths.EPI.epiGS,'PCA5','FC*.mat'));
        [~,result] = system(sentence);
        FCscrubbing = PCA5.FC;
        FCrobust = PCA5.FCrobust;
        VProbust = PCA5.VProbust;
        save(fullfile(paths.EPI.epiGS,'PCA5','FCscrubbing.mat'),'FCscrubbing');
        save(fullfile(paths.EPI.epiGS,'PCA5','FCrobust.mat'),'FCrobust');
        save(fullfile(paths.EPI.epiGS,'PCA5','VProbust.mat'),'VProbust');
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
