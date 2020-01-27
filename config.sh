#!/bin/bash


################################################################################
################################################################################
## GLOBALS & dependencies

# source bash funcs
source ${EXEDIR}/src/func/bash_funcs.sh

## pip install pydicom --user

################################################################################
############################  PATH TO DATA  ###################################

#export path2data="/N/dc2/scratch/aiavenak/testdata"
export path2data="/N/dc2/projects/connectivitypipeline/example_for_andrea/SUBJECTS"

################################################################################
#####################  SET UP DIRECTORY STRUCTURE  #############################

# The following diagrapm is a sample directory tree for a single subject.
# Following that are configs you can use to set your own names if different
# from sample structure.

# SUBJECT1 -- T1 -- DICOMS
#          |
#          -- EPI(#) -- DICOMS (May have multiple EPI scans)
#          |         |
#          |         |               (SPIN-ECHO)       (GRADIENT ECHO)
#          |         -- UNWARP -- SEFM_AP_DICOMS (OR) GREFM_MAG_DICOMS
#          |         |         | 
#          |         |         -- SEFM_PA_DICOMS (OR) GREFM_PHASE_DICOMS
#          |         |
#          |         -- UNWARPED uf*.nii.gz (MELODIC UNWARPED IMAGES)
#          |         |
#          |         |
#          |         -- FEAT -- FEAT_PREP
#          |
#          -- DWI -- DICOMS
#                 |
#                 -- UNWARP -- B0_PA_DCM

export configs_T1="T1"
export configs_epiFolder="EPI"
    export configs_sefmFolder="UNWARP" # Reserved for Field Mapping series
        export configs_APdcm="SEFM_AP_DICOMS" # Spin Echo A-P
        export configs_PAdcm="SEFM_PA_DICOMS" # Spin Echo P-A
        export configs_GREmagdcm="GREFM_MAG_DICOMS" # Gradient echo FM magnitude series
        export configs_GREphasedcm="GREFM_PHASE_DICOMS" # Gradient echo FM phase map series
	# export configs_melodicUnwarpedFolder="UNWARPED"
	# export configs_FEAT="FEAT"	
		# export configs_FEAT_PREP="FEAT_PREP"
export configs_DWI="DWI"
    export configs_unwarpFolder="UNWARP"
        export configs_dcmPA="B0_PA_DCM" #b0 opposite phase encoding

export configs_dcmFolder="DICOMS"
export configs_dcmFiles="IMA.dcm" # Dicom file extension
export configs_niiFiles="nii" # Nifti-1 file extension


################################################################################
################################ TEMPLATES #####################################

export pathFSLstandard="${FSLDIR}/data/standard"
## path to Supplementary Materials (SM)
export pathSM="/N/dc2/projects/connectivitypipeline/example_for_andrea/ConnPipelineSM"
export pathMNItmplates="${pathSM}/MNI_templates"
export pathBrainmaskTemplates="${pathSM}/brainmask_templates"
export pathParcellations="${pathSM}/Parcellations"
export PYpck="${pathSM}/python-pkgs"


################################################################################
################################ PARCELLATIONS #################################

# required parc
export PARC0="CSFvent"
export PARC0dir="${pathMNItmplates}/MNI152_T1_1mm_VentricleMask.nii.gz"
export PARC0pcort=0;
export PARC0pnodal=0;

# required parc
export PARC1="shen_278"
export PARC1dir="shen_MNI152_org"
export PARC1pcort=0;
export PARC1pnodal=1;

# optional
export PARC2="yeo7"
export PARC2dir="yeo7_MNI152"
export PARC2pcort=1;
export PARC2pnodal=0;

# optional
export PARC3="yeo17"
export PARC3dir="yeo17_MNI152"
export PARC3pcort=1;
export PARC3pnodal=0;

# # Schaefer parcellation of yeo17 into 200 nodes
# parcs.plabel(1).name='schaefer200_yeo17';
# parcs.pdir(1).name='Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm';
# parcs.pcort(1).true=1;
# parcs.pnodal(1).true=1;

# # Schaefer parcellation of yeo17 into 300 nodes
# parcs.plabel(2).name='schaefer300_yeo17';
# parcs.pdir(2).name='Schaefer2018_300Parcels_17Networks_order_FSLMNI152_1mm';
# parcs.pcort(2).true=1;
# parcs.pnodal(2).true=1;

export numParcs=2  # CSF doesn't count; numParcs cannot be less than 1. Shen is the defailt parc


################################################################################
############################# T1_PREPARE_A #####################################

export T1_PREPARE_A=true

if $T1_PREPARE_A; then

	export flags_T1_dcm2niix=true;  # dicom to nifti conversion 
		export configs_T1_useCropped=false; # use cropped field-of-view output of dcm2niix
		
	export flags_T1_denoiser=true # denoising
		# export configs_T1_denoised="T1_denoised_SUSAN"  ## this should eventually be an input param SUSAN vs ANTS
	
	export flags_T1_anat=true # run FSL_anat
		export configs_T1_bias=1; # 0 = no; 1 = weak; 2 = strong
		export configs_T1_crop=0; # 0 = no; 1 = yes (lots already done by dcm2niix)

	export flags_T1_bet=true; # brain extraction and mask generation (only needed for double BET)
		export configs_antsTemplate="MICCAI"  # options are: MICCAI, NKI or bet
		export configs_T1_A_betF="0.35" # this are brain extraction parameters with FSL bet
		export configs_T1_A_betG="0.15"  # see fsl bet help page for more details
		# ANTS does not require bet inputs
	 
	export flags_T1_re_extract=true; # brain extraction with mask

fi 

# Set denoising option
export flag_ANTS=true # other option available is FSL's SUSAN, set flag_ANTS=false to use SUSAN instead 
if ${flag_ANTS}; then 
	export configs_T1_denoised="T1_denoised_ANTS" 
else
	export configs_T1_denoised="T1_denoised_SUSAN"
fi


################################################################################
############################# T1_PREPARE_B #####################################

export T1_PREPARE_B=true

if $T1_PREPARE_B; then

	# registration flags
	export flags_T1_reg2MNI=true
		export configs_T1_useExistingMats=false
		export configs_T1_useMNIbrain=false
		export configs_T1_fnirtSubSamp="4,4,2,1"
	# segmentation flags
	export flags_T1_seg=true		
		export configs_T1_segfastH="0.25"
		export configs_T1_masklowthr=1
		export configs_T1_flirtdof6cost="mutualinfo"
	# parcellation flags
	export flags_T1_parc=true
		export configs_T1_numDilReMask=3
		export configs_T1_addsubcort=true # ad FSL subcortical to cortial parcellations ONLY

	export path2MNIref="${pathFSLstandard}/MNI152_T1_1mm.nii.gz"

fi 


################################################################################
############################# fMRI_A #####################################

export fMRI_A=true

if $fMRI_A; then

	# # set number of EPI sessions/scans
	export configs_EPI_epiMin=1; # minimum scan index
	export configs_EPI_epiMax=4; # maximum scan index

	export flags_EPI_dcm2niix=true; # dicom import

	export flags_EPI_ReadHeaders=true; # obtain pertinent scan information
		export flags_EPI_UseJson=true; # obtain pertinent scan information through json files generated by dcm2niix
		export scanner_param_TR="RepetitionTime"  # "RepetitionTime" for Siemens; "tr" for GE
		export scanner_param_TE="EchoTime"  # "EchoTime" for Siemens; "te" for GE
		export scanner_param_FlipAngle="FlipAngle"  # "FlipAngle" for Siemens; "flip_angle" for GE
		export scanner_param_EffectiveEchoSpacing="EffectiveEchoSpacing"  # "EffectiveEchoSpacing" for Siemens; "effective_echo_spacing" for GE
		export scanner_param_BandwidthPerPixelPhaseEncode="BandwidthPerPixelPhaseEncode"  # "BandwidthPerPixelPhaseEncode" for Siemens; unknown for GE
		export scanner_param_slice_fractimes="SliceTiming"  # "SliceTiming" for Siemens; "slice_timing" for GE

	export flags_EPI_SpinEchoUnwarp=true # Requires UNWARP directory and approporiate dicoms.
	# # SPIN ECHO PAIRS (A-P, P-A) Acquistion on the Prisma
		export configs_EPI_SEnumMaps=3; # Fallback Number of PAIRS of AP and PA field maps.
	# # Defaults to reading *.dcm/ima files in SE AP/PA folders
	# # topup (see www.mccauslanddenter.sc.edu/cml/tools/advanced-dti - Chris Rorden's description
	# # readOutTime=echoSpacing*((matrixlines4phase*partialFourier/accelrerationFactor)-1)
	# # readOutTime now calculated from image data.
		export flags_EPI_RunTopup=true # 1=Run topup (1st pass), 0=Do not rerun if previously completed.       
	# # Gradient recalled echo Field Map Acquisition
	# configs.EPI.GREmagdcm='GREFM_MAG_DICOMS'; # MAGNITUDE Series
	# configs.EPI.GREphasedcm='GREFM_PHASE_DICOMS'; # PHASE Series
	# configs.EPI.GREbetf=0.5; # GRE-specific bet values. Do not change
	# configs.EPI.GREbetg=0;   # GRE-specific bet input. Change if needed 
	# configs.EPI.GREdespike=1; # Perform FM despiking
	# configs.EPI.GREsmooth=3; # GRE phase map smoothing (Gaussian sigma, mm)
	# configs.EPI.EPIdwell=0.000308; # Dwell time (sec) for the EPI to be unwarped 
	export flags_EPI_SliceTimingCorr=true
		#export flags_EPI_UseUnwarped=true # Use unwarped EPI if both warped and unwarped are available.
		export configs_EPI_UseTcustom=1;# 1: use header-extracted times (suggested)

	export flags_EPI_MotionCorr=true

	export flags_EPI_RegT1=true;
		export configs_EPI_epibetF=0.3000;
		export configs_EPI_minVoxelsClust=8; # originally hardwired to 8

	export flags_EPI_RegOthers=true;
		export configs_EPI_GMprobthr=0.2; # Threshold the GM probability image
										# change from 0.25 to 0.2 or 0.15

	export flags_EPI_IntNorm4D=false; # Intensity normalization to global 4D mean of 1000

	export flags_EPI_AROMA=false; # ICA-based denoising; WARNING: This will smooth your data.
		export ICA_AROMA="${PYpck}/ICA-AROMA"
		if [[ -e "${pathFSLstandard}/MNI152_T1_2mm_brain.nii.gz" ]]; then
			fileMNI2mm="${pathFSLstandard}/MNI152_T1_2mm_brain.nii.gz"
		else
			fileMNI2mm="${pathMNItmplates}/MNI152_T1_2mm_brain.nii.gz"
		fi


	export flags_EPI_DemeanDetrend=false;

	export flags_EPI_MotionRegressors=false
		export configs_EPI_scrubtime=15
	# GS is a subflag of Motion Regressors. 
	# If equals 1 global signal regression is done, if 0 it is not done.
		export flags_EPI_GS=1
		export configs_EPI_FDth='0.20';

fi

################################################################################
############################# DWI processing ###################################

export DWI_A=false

if $DWI_A; then

	export flags_DWI_dcm2niix=false # dicom to nifti coversion
		export configs_DWI_readout=[] # if empty get from dicom; else specify value
	export flags_DWI_topup=false # FSL topup destortion field estimation
		export configs_DWI_b0cut=1 # maximum B-value to be considered B0
	export flags_DWI_eddy=false # FSL EDDY distortion correction
		export configs_DWI_EDDYf='0.3' # fsl bet threshold for b0 brain mask used by EDDY
		export configs_DWI_repolON=true # use eddy_repol to interpolate missing/outlier data
	export flags_DWI_DTIfit=true  # Tensor estimation and generation of scalar maps
		export configs_DWI_DTIfitf='0.4' # brain extraction (FSL bet -f) parameter 

fi 

