
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


echo "# =========================================================="
echo "# 4. Normalization to 4D mean 1000."
echo "# =========================================================="

if [[ ! -e "${EPIpath}/3_epi.nii.gz" ]]; then  
    log "WARNING File ${EPIpath}/3_epi.nii.gz does not exist. Skipping further analysis"
    exit 1 
else
    fileIn="${EPIpath}/3_epi.nii.gz"
    fileOut="${EPIpath}/4_epi.nii.gz"

    cmd="fslmaths ${fileIn} -ing 1000 ${fileOut}"
    log $cmd
    eval $cmd 
fi