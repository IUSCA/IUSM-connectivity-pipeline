               
#!/bin/bash
#
# Script: fMRI_A adaptaion from Matlab script 
#

###############################################################################
#
# Environment set up
#
###############################################################################

shopt -s nullglob # No-match globbing expands to null

source ${EXEDIR}/src/func/bash_funcs.sh

############################################################################### 

            
echo "# ==========================================="
echo "# 2. Motion Correction"
echo "# ==========================================="

echo "**** ${EPIpath}/1_epi.nii.gz"

if [[ ! -e "${EPIpath}/1_epi.nii.gz" ]]; then  
    log "-No slice time corrected 1_epi output found"
    log "  Defaulting to 0_epi data."
    if [[ ! -e "${EPIpath}/0_epi_unwarped.nii.gz" ]]; then 
        log " -Unwarped 0_epi volume does not exist"
        fileIn="${EPIpath}/0_epi.nii.gz"
        if [[ -e "${fileIn}" ]]; then 
            log " -Will use 0_epi from dicom conversion."
        else
            log "WARNING No 0_epi inputs found... Exiting"
            exit 1
        fi
    else
        fileIn="${EPIpath}/0_epi_unwarped.nii.gz"
        log " -Will use 0_epi_unwarped.nii.gz"
    fi 
else
    fileIn="${EPIpath}/1_epi.nii.gz"
    log " -Will use the slice time corrected 1_epi.nii.gz as input" 
fi 

# Compute motion outliers
cmd="${EXEDIR}/src/scripts/get_motion_outliers.sh ${EPIpath} ${fileIn}"
log $cmd
eval $cmd


fileOut="${EPIpath}/2_epi"
cmd="fslval ${fileIn} dim4"
log $cmd 
nvols=`$cmd`                
#nvols=`echo $out | awk -F' ' '{ print $2}'`
echo "export nvols=${nvols}" >> ${EPIpath}/0_param_dcm_hdr.sh
echo "Number of volumes in 1_epi_brain: ${nvols} "

cmd="mcflirt -in ${fileIn} -out ${fileOut} -plots -meanvol"
log $cmd 
eval $cmd 

cmd="mv ${EPIpath}/2_epi.par ${EPIpath}/motion.txt"
log $cmd 
eval $cmd 