
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

function read_data() {
fileIN="$1" EPIpath="$2" ${python3_7} - <<END
import os
import numpy as np
import nibabel as nib


fileIN=os.environ['fileIN']
EPIpath=os.environ['EPIpath']

# load data and masks
resting = nib.load(fileIN)
resting_vol = resting.get_data()
[sizeX,sizeY,sizeZ,numTimePoints] = resting_vol.shape
print(sizeX,sizeY,sizeZ,numTimePoints)

fname = ''.join(EPIpath,'/rT1_CSFvent_mask_eroded.nii.gz'])
volCSFvent = nib.load(fname)
volCSFvent_vol = volCSFvent.get_data()

END
}

##############################################################################

## PHYSIOLOGICAL REGRESSORS
log "============== 5.1 PHYSIOLOGICAL REGRESSORS =================="

if ${flags_PhysiolReg_aCompCorr}; then  
    log "- -----------------aCompCor---------------------"      
elif ${flags_PhysiolReg_WM_CSF}; then
    log "- -----------------Mean CSF and WM Regression-----------------"
fi 

if ${flags_NuisanceReg_AROMA}; then   

    fileIN="${EPIpath}/AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz"
    if  [[ -e ${fileIN} ]]; then
        log "PhysiolReg - Using AROMA output data"
        if ${flags_PhysiolReg_aCompCorr}; then  
            PhReg_path="${EPIpath}/AROMA/aCompCorr"    
        elif ${flags_PhysiolReg_WM_CSF}; then
            PhReg_path="${EPIpath}/AROMA/PhysReg"
        fi          
    else
        log "ERROR ${fileIN} not found. Connot perform physiological regressors analysis"
    fi 

elif ${flags_NuisanceReg_HeadParam}; then 

    fileIN="${EPIpath}/4_epi.nii.gz"
    HMPpath="${EPIpath}/HMPreg"
    if  [[ -e ${fileIN} ]] && [[ -d ${HMPpath} ]]; then
        if ${flags_PhysiolReg_aCompCorr}; then  
            log "PhysiolReg - Combining aCompCorr with HMP regressors"
            PhReg_path="${EPIpath}/HMPreg/aCompCorr"    
        elif ${flags_PhysiolReg_WM_CSF}; then
            log "PhysiolReg - Combining Mean CSF & WM signal with HMP regressors"
            PhReg_path="${EPIpath}/HMPreg/PhysReg"
        fi          
    else
        log "ERROR ${fileIN} and or ${HMPpath} not found. Connot perform physiological regressors analysis"
    fi 
fi


if [[ ! -d ${PhReg_path} ]]; then
    cmd="mkdir ${PhReg_path}"
    log $cmd
    eval $cmd 
fi

# read in data and masks 
read_data ${fileIN}
