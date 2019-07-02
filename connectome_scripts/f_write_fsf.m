function [paths, configs, params] = f_write_fsf(configs, params, paths)
%write_fsf Writes a FEAT structure file (.fsf) for running Melodic single
%session.
%   fsf file is a structure set up file for running FEAT, including
%   Melodic, usually it can be saved within the Melodic GUI interface, or
%   we can write it ourselves, for more specific setup.

% Contributors:
%   Zikai Lin, Indiana University School of Medicine
%
% Updates:
%   03/26/2019 Initial update. -ZL


% Open an fsf file in subject's EPI directory
fsf_path = fullfile( paths.EPI.feat_prep_dir, sprintf('%d_%s_SingleMELODIC.fsf',configs.EPI.currEPI, configs.EPI.subject));
fid = fopen(fsf_path,'wt');
configs.EPI.subjectFSF = fsf_path;
configs.EPI.feat_subject_out = fullfile(paths.EPI.feat_out, configs.EPI.subject);
% FEAT version number
fsf_str = sprintf([ ...
       '\n# FEAT version number' ...
       '\nset fmri(version) %s\n'], configs.EPI.featVersion);
fprintf(fid,'%s', fsf_str);

% Melodic or not, here set to 1
fsf_str = sprintf([ ...
       '\n# Are we in MELODIC?' ...
       '\nset fmri(inmelodic) 1\n']);
fprintf(fid,'%s', fsf_str);

% Analysis Level
fsf_str = sprintf([ ...
       '\n# Analysis level\n# 1 : First-level analysis\n# 2 : Higher-level analysis' ...
       '\nset fmri(level) 1\n']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Which stages to run, here we run 
fsf_str = sprintf([ ...
       '# Which stages to run\n# 0 : No first-level analysis (registration and/or group stats only)',...
       '\n# 7 : Full first-level analysis',...
       '\n# 1 : Pre-processing',...
       '\n# 2 : Statistics',...
       '\nset fmri(analysis) 7']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Use relative filenames or not, here we set to 0
fsf_str = sprintf([ ...
       '\n# Use relative filenames',...
       '\nset fmri(relative_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Use balloon help or not, here we set to 1
fsf_str = sprintf([ ...
       '\n# Balloon help',...
       '\nset fmri(help_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Use Featwatcher or not, this is up to user
% Featwatcher is a simple GUI that presents the progress of a FEAT analysis.
% When you start FEAT running with the Go button, it automatically brings up the Featwatcher GUI.
% At the top is the name of the FEAT output directory; 
% under that is a small panel which shows some running statistics on the 
% currently running FEAT sub-process - %CPU usage, total CPU time so far 
% and memory usage. In the large panel is shown the report.log file which 
% lists all commands run by FEAT, along with their output.
fsf_str = sprintf([ ...
       '\n# Run Featwatcher',...
       '\nset fmri(featwatcher_yn) %d'], configs.EPI.watcher);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Cleanup first-level standard-space images, here we set to 0
fsf_str = sprintf([ ...
       '\n# Cleanup first-level standard-space images',...
       '\nset fmri(sscleanup_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Cleanup first-level standard-space images, here we set to 0
fsf_str = sprintf([ ...
       '\n# Output directory',...
       '\nset fmri(outputdir) %s'], fullfile(paths.EPI.feat_out, configs.EPI.subject));
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% TR, can be read from dicom header or json file
fsf_str = sprintf([ ...
       '\n# TR(s)',...
       '\nset fmri(tr) %f'], params.TR);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Total volumes, can be read from dicom header or json file
fsf_str = sprintf([ ...
       '\n# Total volumes',...
       '\nset fmri(npts) %d'], params.nvols);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Delete volumes, set by user in batch_set_up.m
fsf_str = sprintf([ ...
       '\n# Delete volumes',...
       '\nset fmri(ndelete) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Perfusion tag/control order, set to 1 by default
fsf_str = sprintf([ ...
       '\n# Perfusion tag/control order',...
       '\nset fmri(tagfirst) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Number of first-level analyses

% Use First-level analysis for analysing each session's data - 
% i.e. the time-series analysis of the raw 4D FMRI data.

% Use Higher-level analysis for combining first-level analyses.
% Here we set to 1 because we run each subject at a time
fsf_str = sprintf([ ...
       '\n# Number of first-level analyses',...
       '\nset fmri(multiple) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Higher level input type
fsf_str = sprintf([ ...
       '\n# Higher-level input type',...
       '\n# 1 : Inputs are lower-level FEAT directories',...
       '\n# 2 : Inputs are cope images from FEAT directories',...
       '\nset fmri(inputtype) 2']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Pre-stat processing?
fsf_str = sprintf([ ...
       '\n# Carry out pre-stats processing?',...
       '\nset fmri(filtering_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Brain background threshold %
fsf_str = sprintf([ ...
       '\n# Brain background threshold, in percentage',...
       '\nset fmri(brain_thresh) %f'], configs.EPI.brainThres);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Critical z for design efficiency calculation, default set to be 5.3
fsf_str = sprintf([ ...
       '\n# Critical z for design efficiency calculation',...
       '\nset fmri(critical_z) 5.3']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Noise level, don't need to be changed
fsf_str = sprintf([ ...
       '\n# Noise level',...
       '\nset fmri(noise) 0.66']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Noise AR
fsf_str = sprintf([ ...
       '\n# Noise AR(1)',...
       '\nset fmri(noisear) 0.34']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Motion Correction, set to be 1, see 
% <https://doi.org/10.1016/j.neuroimage.2015.02.064>
fsf_str = sprintf([ ...
       '\n# Motion correction',...
       '\n# 0 : None',...
       '\n# 1 : MCFLIRT'...
       '\nset fmri(mc) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Spin History (currently obsolete) just let it set to be default
fsf_str = sprintf([ ...
       '\n# Spin-history (currently obsolete)',...
       '\nset fmri(sh_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% B0 fieldmap unwarping
fsf_str = sprintf([ ...
       '\n# B0 fieldmap unwarping?',...
       '\nset fmri(regunwarp_yn) %d'], configs.EPI.B0Unwarp);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% GDC Test, set to be None by default
fsf_str = sprintf([ ...
       '\n# GDC Test',...
       '\nset fmri(gdc) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% EPI dwell time (ms), can be set in batch_set_up.m, we are not using
% unwarping in melodic prior to ICA-AROMA, so the EPIdwell time doesn't
% matter.
fsf_str = sprintf([ ...
       '\n# EPI dwell time (ms)',...
       '\nset fmri(dwell) %f'], configs.EPI.EPIdwell*1000);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% EPI TE (ms), can be read from dicom headers or json files
fsf_str = sprintf([ ...
       '\n# EPI TE (ms)',...
       '\nset fmri(te) %f'], params.TE*1000);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% % Signal loss threshold, leave this as default
fsf_str = sprintf([ ...
       '\n# Percentage of Signal loss threshold',...
       '\nset fmri(signallossthresh) 10']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Unwarp direction
fsf_str = sprintf([ ...
       '\n# Unwarp direction',...
       '\nset fmri(unwarp_dir) y-']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Slice timing correction, can be done later in the pipeline, so leave it 0
% here
fsf_str = sprintf([ ...
       '\n# Slice timing correction',...
       '\n# 0 : None',...
       '\n# 1 : Regular up (0,1,2,3,...)',...
       '\n# 2 : Regular down',...
       '\n# 3 : Use slice order file',...
       '\n# 4 : Use slice timings file',...
       '\n# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )',...
       '\nset fmri(st) %d'], configs.EPI.melodic_st);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% # Slice timings file
fsf_str = sprintf([ ...
       '\n# Unwarp direction# Slice timings file',...
       '\nset fmri(st_file) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% BET brain extraction, we need it so set to 1
fsf_str = sprintf([ ...
       '\n# BET brain extraction',...
       '\nset fmri(bet_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


%  Spatial smoothing FWHM (mm), if this has been done, it will skip the
%  spatial smoothing later in the pipeline.
fsf_str = sprintf([ ...
       '\n# Spatial smoothing FWHM (mm)',...
       '\nset fmri(smooth) %f'], configs.EPI.pre_fwhm);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% # Intensity normalization, hardwired to 1 (yes)
fsf_str = sprintf([ ...
       '\n# Intensity normalization',...
       '\nset fmri(norm_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Perfusion subtraction
fsf_str = sprintf([ ...
       '\n# Perfusion subtraction',...
       '\nset fmri(perfsub_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Highpass temporal filtering
fsf_str = sprintf([ ...
       '\n# Highpass temporal filtering',...
       '\nset fmri(temphp_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Lowpass temporal filtering, it is not useful so this is turned of by
% Melodic. See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide#First-level_or_Higher-level_Analysis.3F
fsf_str = sprintf([ ...
       '\n# Lowpass temporal filtering',...
       '\nset fmri(templp_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% MELODIC ICA data exploration, default turn off
fsf_str = sprintf([ ...
       '\n# MELODIC ICA data exploration',...
       '\nset fmri(melodic_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Carry out main stat? yes or no.
fsf_str = sprintf([ ...
       '\n# Carry out main stats?',...
       '\nset fmri(stats_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Carry out prewhitening? According to Melodic user manual, prewhitening
% should be turned off if TR > 30s or number of time points(ntp) < 50.
if (params.TR > 30 || params.nvols < 50)
    fsf_str = sprintf([ ...
       '\n# Carry out prewhitening?',...
       '\nset fmri(prewhiten_yn) 0']);
    fprintf(fid,'%s', fsf_str);
else
    fsf_str = sprintf([ ...
       '\n# Carry out prewhitening?',...
       '\nset fmri(prewhiten_yn) 1']);
    fprintf(fid,'%s', fsf_str);
end
fprintf(fid,'\n');


% Add motion parameters to model
fsf_str = sprintf([ ...
       '\n# Add motion parameters to model',...
       '\n# 0 : No',...
       '\n# 1 : Yes',...
       '\n set fmri(motionevs) 0',...
       '\n set fmri(motionevsbeta) ""',...
       '\n set fmri(scriptevsbeta) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

%  Robust outlier detection in FLAME? Set to 0 by default
fsf_str = sprintf([ ...
       '\n# Robust outlier detection in FLAME?',...
       '\nset fmri(robust_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Higher-level modeling, irrelevant to MELODIC, set it as default
fsf_str = sprintf([ ...
       '\n# Higher-level modeling',...
       '\n# 3 : Fixed effects',...
       '\n# 0 : Mixed Effects: Simple OLS',...
       '\n# 2 : Mixed Effects: FLAME 1',...
       '\n# 1 : Mixed Effects: FLAME 1+2',...
       '\nset fmri(mixed_yn) 2']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Higher-level permutations, set it to 5000 as default
fsf_str = sprintf([ ...
       '\n# Higher-level permutations',...
       '\nset fmri(randomisePermutations) 5000']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Number of EVs, irrelevant to MELODIC, set it as default
fsf_str = sprintf([ ...
       '\n# Number of EVs',...
       '\nset fmri(evs_orig) 1',...
       '\nset fmri(evs_real) 1',...
       '\nset fmri(evs_vox) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Number of contrasts, irrelevant to MELODIC, set it as default
fsf_str = sprintf([ ...
       '\n# Number of contrasts',... 
       '\nset fmri(ncon_orig) 1',...
       '\nset fmri(ncon_real) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Number of F-test, irrelevant to MELODIC, set it as default
fsf_str = sprintf([ ...
      '\n# Number of F-test',...     
       '\nset fmri(nftests_orig) 0',...
       '\nset fmri(nftests_real) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Add constant column to design matrix? (obsolete)
fsf_str = sprintf([ ...
       '\n# Add constant column to design matrix? (obsolete)',...   
       '\nset fmri(constcol) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Carry out post-stats steps?
fsf_str = sprintf([ ...
       '\n# Carry out post-stats steps?',...   
       '\nset fmri(poststats_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Pre-threshold masking?
fsf_str = sprintf([ ...
       '\n# Pre-threshold masking?',...   
       '\nset fmri(threshmask) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Thresholding
fsf_str = sprintf([ ...
       '\n#Thresholding\n# 0 : None\n# 1 : Uncorrected\n# 2 : Voxel\n# 3 : Cluster',...   
       '\nset fmri(thresh) 3']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% P threshold
fsf_str = sprintf([ ...
       '\n# P threshold',...   
       '\nset fmri(prob_thresh) 0.05']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Z threshold, set as default.
fsf_str = sprintf([ ...
       '\n# Z threshold',...   
       '\nset fmri(z_thresh) 3.1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Z min in colour rendering
fsf_str = sprintf([ ...
       '\n# 0 : Use actual Z min/max \n# 1 : Use preset Z min/max',...   
       '\nset fmri(zdisplay) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Z min in colour rendering
fsf_str = sprintf([ ...
       '\n# Z min in colour rendering',...   
       '\nset fmri(zmin) 2']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Z max in colour rendering
fsf_str = sprintf([ ...
       '\n# Z max in colour rendering',...   
       '\nset fmri(zmax) 8']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Colour rendering type
fsf_str = sprintf([ ...
       '\n# Colour rendering type\n# 0 : Solid blobs\n# 1 : Transparent blobs',...   
       '\nset fmri(rendertype) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Background image for higher-level stats overlays
fsf_str = sprintf([ ...
       '\n# Background image for higher-level stats overlays\n# 1 : Mean highres\n# 2 : First highres\n# 3 : Mean functional\n# 4 : First functional\n# 5 : Standard space template',...   
       '\nset fmri(bgimage) %d'], configs.EPI.bgimage);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Create time series plots
fsf_str = sprintf([ ...
       '\n# Create time series plots',...   
       '\nset fmri(tsplot_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Registration to initial structural
fsf_str = sprintf([ ...
       '\n# Registration to initial structural',...   
       '\nset fmri(reginitial_highres_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Search space for registration to initial structural
fsf_str = sprintf([ ...
       '\n#Search space for registration to initial structural\n# 0   : No search\n# 90  : Normal search\n# 180 : Full search',...   
       '\nset fmri(reginitial_highres_search) 90']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Degrees of Freedom for registration to initial structural
fsf_str = sprintf([ ...
       '\n# Degrees of Freedom for registration to initial structural',...   
       '\nset fmri(reginitial_highres_dof) 3']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Registration to main structural, set as yes by default
fsf_str = sprintf([ ...
       '\n# Registration to main structural',...   
       '\nset fmri(reghighres_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Search space for registration to main structural, set as normal search.
fsf_str = sprintf([ ...
       '\n# Search space for registration to main structural\n# 0   : No search\n# 90  : Normal search\n# 180 : Full search',...   
       '\nset fmri(reghighres_search) %d'], configs.EPI.reghighres_search);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Degrees of Freedom for registration to main structural
fsf_str = sprintf([ ...
       '\n# Degrees of Freedom for registration to main structural',...   
       '\nset fmri(reghighres_dof) BBR']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Registration to standard image?
fsf_str = sprintf([ ...
       '\n# Registration to standard image?',...   
       '\nset fmri(regstandard_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Use alternate reference images?
fsf_str = sprintf([ ...
       '\n# Use alternate reference images?',...   
       '\nset fmri(alternateReference_yn) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Standard image
fsf_str = sprintf([ ...
       '\n# Standard image',...   
       '\nset fmri(regstandard) %s'], paths.stdImg);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

%  Search space for registration to standard space
fsf_str = sprintf([ ...
       '\n# Search space for registration to standard space\n# 0   : No search\n# 90  : Normal search\n# 180 : Full search',...   
       '\nset fmri(regstandard_search) %d'], configs.EPI.regstandard_search);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Degrees of Freedom for registration to standard space
fsf_str = sprintf([ ...
       '\n# Degrees of Freedom for registration to standard space',...   
       '\nset fmri(regstandard_dof) %d'], configs.EPI.regstandard_dof);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Do nonlinear registration from structural to standard space?
fsf_str = sprintf([ ...
       '\n# Do nonlinear registration from structural to standard space?',...   
       '\nset fmri(regstandard_nonlinear_yn) %d'], configs.EPI.regstandard_nonlinear_yn);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Control nonlinear warp field resolution (mm)
fsf_str = sprintf([ ...
       '\n# Control nonlinear warp field resolution',...   
       '\nset fmri(regstandard_nonlinear_warpres) %d'], configs.EPI.regstandard_nonlinear_warpres);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% High pass filter cutoff (s)
fsf_str = sprintf([ ...
       '\n# High pass filter cutoff',...   
       '\nset fmri(paradigm_hp) 100'], configs.EPI.paradigm_hp);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Total voxels
% ?Total voxels? is used purely for estimating the run-time when submitting to a cluster -
% you can safely set it to an arbitrary default value in a template design.fsf.
fsf_str = sprintf([ ...
       '\n# Total voxels',...   
       '\nset fmri(totalVoxels) 202954752']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Number of lower-level copes feeding into higher-level analysis
fsf_str = sprintf([ ...
       '\n# Number of lower-level copes feeding into higher-level analysis',...   
       '\nset fmri(ncopeinputs) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% 4D AVW data or FEAT directory (1)
fsf_str = sprintf([ ...
       '\n# 4D AVW data or FEAT directory (1)',...   
       '\nset feat_files(1) %s'], paths.EPI.feat_4Dnifti);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Add confound EVs text file
fsf_str = sprintf([ ...
       '\n# Add confound EVs text file',...   
       '\nset fmri(confoundevs) 0']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Subjects structural image for analysis 
fsf_str = sprintf([ ...
       '\n# Subjects structural image for analysis 1',...   
       '\nset highres_files(1) %s'], paths.EPI.structImage);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Resampling resolution
fsf_str = sprintf([ ...
       '\n# Resampling resolution',...   
       '\nset fmri(regstandard_res) %d'], configs.EPI.regstandard_res);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Variance-normalise timecourses
fsf_str = sprintf([ ...
       '\n# Variance-normalise timecourses',...   
       '\nset fmri(varnorm) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Automatic dimensionality estimation
fsf_str = sprintf([ ...
       '\n# Automatic dimensionality estimation',...   
       '\nset fmri(dim_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Output components
fsf_str = sprintf([ ...
       '\n# Output components',...   
       '\nset fmri(dim) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Types of ICA, sigle session melodic as default
fsf_str = sprintf([ ...
       '\n# Types of ICA \n# 1 : Single-session ICA\n# 2 : Multi-session temporal concatenation\n# 3 : Multi-session tensor TICA',...   
       '\nset fmri(icaopt) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Threshold IC maps
fsf_str = sprintf([ ...
       '\n# Threshold IC maps',...   
       '\nset fmri(thresh_yn) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Mixture model threshold
fsf_str = sprintf([ ...
       '\n# Mixture model threshold',...   
       '\nset fmri(mmthresh) %f'], configs.EPI.mmthresh);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Output full stats folder
fsf_str = sprintf([ ...
       '\n# Output full stats folder',... 
       '\nset fmri(ostats) 1']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Timeseries and subject models, use default setting
fsf_str = sprintf('\n# Timeseries and subject models\nset fmri(ts_model_mat) ""\nset fmri(ts_model_con) ""\nset fmri(subject_model_mat) ""\nset fmri(subject_model_con) ""');
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid, "%s\n", "##########################################################");
fprintf(fid, "# Now options that don't appear in the GUI\n");

%% Options that dont appear in the GUI

% Alternative (to BETting) mask image
fsf_str = sprintf([ ...
       '\n# Alternative (to BETting) mask image',...   
       '\nset fmri(alternative_mask) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Initial structural space registration initialisation transform
fsf_str = sprintf([ ...
       '\n# Initial structural space registration initialisation transform',...   
       '\nset fmri(init_initial_highres) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');


% Structural space registration initialisation transform
fsf_str = sprintf([ ...
       '\n# Structural space registration initialisation transform',...   
       '\nset fmri(init_highres) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% Standard space registration initialisation transform
fsf_str = sprintf([ ...
       '\n# Standard space registration initialisation transform',...   
       '\nset fmri(init_standard) ""']);
fprintf(fid,'%s', fsf_str);
fprintf(fid,'\n');

% For full FEAT analysis: overwrite existing .feat output dir?
fsf_str = sprintf([ ...
       '\n# For full FEAT analysis: overwrite existing .feat output dir?',...   
       '\nset fmri(overwrite_yn) 0']);
fprintf(fid,'%s', fsf_str);

fclose(fid);
end

