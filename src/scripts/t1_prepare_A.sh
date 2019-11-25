
#!/bin/bash
#
# Script: T1_PREPARE_A adaptaion from Matlab script 
#

###############################################################################
#
# Environment set up
#
###############################################################################

shopt -s nullglob # No-match globbing expands to null

source ${EXEDIR}/src/func/bash_funcs.sh



if [[ -d "$T1path/DICOMS" ]]; then
	if [[ "$(ls -A $T1path/DICOMS)" ]]; then 

        log "T1_PREPARE_A"

		##### DICOM 2 nifti ######

		if ${flags_T1_dcm2niix}; then
			# if .IMA or .IMA.dcm or .dcm files exis inside T1/DICOMS
			niifiles=`find $T1path -maxdepth 1 -name "*.nii*" | wc -l`
			if [[ $niifiles != 0 ]]; then 
				log "rm $T1path/*.nii"
				rm $T1path/*.nii*
			fi
			
			cmd="dcm2niix -f T1 -o $T1path -v y -x y $T1path/DICOMS > $T1path/dcm2niix.log"
			log $cmd
			eval $cmd 
			if [[ $? != 0 ]];  then
				echoerr "dcm2nii exit status not zero. exiting."
				exit 1
			else
				cmd="mv $T1path/T1.nii $T1path/T1_orig.nii"
				log $cmd
				eval $cmd 
				if ${configs_T1_useCropped}; then
					gzip -f $T1path/T1_Crop_1.nii
					mv $T1path/T1_Crop_1.nii.gz $T1path/T1.nii.gz
				else
					gzip -f $T1path/T1_orig.nii
					mv $T1path/T1_orig.nii.gz $T1path/T1.nii.gz
				fi
			fi			
		fi

		##### T1 denoiser ######
		
		if ${flags_T1_denoiser}; then
			#singularity pull --name ants.simg docker://brainlife/ants 
			#singularity shell ants.simg 
			if ${flag_ANTS}; then 
				cmd="DenoiseImage -v -d 3 -n Gaussian -p 1 -r 1 -i $T1path/T1.nii.gz -o $T1path/${configs_T1_denoised}.nii.gz"
			else
				cmd="susan $T1path/T1 56.5007996 3 3 1 0 $T1path/${configs_T1_denoised}"
			fi			
			log $cmd
			eval $cmd 
		fi 

		##### FSL ANAT ######
		if ${flags_T1_anat}; then

			if [ ${configs_T1_bias} -eq "0" ]; then
				T1bias="--nobias"
			elif [ ${configs_T1_bias} -eq "1" ]; then
				T1bias="--weakbias"
			else 
				T1bias=" "
			fi


			if [ ${configs_T1_crop} -eq "0" ]; then
				T1crop="--nocrop"
			else 
				T1crop=" "
			fi


			T1args="${T1bias} ${T1crop}"

			if [[ -d "${T1path}/${configs_T1_denoised}.anat" ]]; then
				cmd="rm -fr ${T1path}/${configs_T1_denoised}.anat"
				log $cmd
				eval $cmd 
			fi

			if [[ -e "$T1path/${configs_T1_denoised}.nii.gz" ]]; then
				cmd="fsl_anat --noreg --nononlinreg --noseg ${T1args} -i ${T1path}/${configs_T1_denoised}.nii.gz"
				log ${cmd}
				eval ${cmd}
			fi 

			if [[ -e "${T1path}/${configs_T1_denoised}.anat/T1_biascorr.nii.gz" ]]; then
				cmd="cp ${T1path}/${configs_T1_denoised}.anat/T1_biascorr.nii.gz ${T1path}/T1_fov_denoised.nii.gz"
				log $cmd
				eval $cmd 
				gunzip ${T1path}/T1_fov_denoised.nii.gz
				if [[ $? != 0 ]]; then
					echo 'FSL_ANAT NOT COMPLETED'
				fi
			fi 
		fi 

		##### T1 BET ######
		if ${flags_T1_bet}; then
			if [[ -e "$T1path/T1_fov_denoised.nii" ]]; then
				cmd="bet ${T1path}/T1_fov_denoised.nii ${T1path}/T1_brain.nii.gz \
				-m -f ${configs_T1_A_betF} -g ${configs_T1_A_betG}"
				log $cmd
				eval $cmd 

				cmd="fslmaths ${T1path}/T1_brain_mask.nii.gz -fillh ${T1path}/T1_brain_mask_filled.nii.gz"
				log $cmd
				eval $cmd				
			fi
		fi

		##### T1 Brain Re-Extract ######
		if ${flags_T1_re_extract}; then
			if [[ -e "$T1path/T1_fov_denoised.nii" ]] && [[ -e "$T1path/T1_brain_mask_filled.nii.gz" ]]; then
				cmd="fslmaths ${T1path}/T1_fov_denoised.nii -mul ${T1path}/T1_brain_mask_filled.nii.gz ${T1path}/T1_brain.nii.gz"
				log $cmd
				eval $cmd
				cmd="fslmaths ${T1path}/T1_brain_mask.nii.gz -fillh ${T1path}/T1_brain_mask_filled.nii.gz"
				log $cmd 
				eval $cmd 
			fi
		fi

	else
		msg="WARNING T1 directory is empty; skipping subject $SUBJ"
		log $msg
		exit 1
	fi
else
	msg="TWARNING 1 directory doesn't exist; skipping T1_prepare_A for subject $SUBJ"
	log $msg
	exit 1
fi 