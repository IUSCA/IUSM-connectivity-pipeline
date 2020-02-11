
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

function apply_reg() {
EPIpath="$1" nuisanceReg="$2" numPC="$3" numReg="$4" numGS="$5" physReg="$6" ${python3_7} - <<END
import os
import numpy as np
import nibabel as nib

print("inside Python script")


EPIpath=os.environ['EPIpath']
nuisanceReg=os.environ['nuisanceReg']
print("nuisanceReg",nuisanceReg)
numPC=(os.environ['numPC'])
print("numPC",numPC)
numReg=(os.environ['numReg'])
print("numReg",numReg)
numGS=(os.environ['numGS'])
print("numGS",numGS)
physReg=os.environ['physReg']
print("physReg",physReg)

print("REGRESSORS -- Creating regressor matrix with the follwing:")

if nuisanceReg == "AROMA":
    print("Applying AROMA regressors")
    resting_file = ''.join([EPIpath,'AROMA/AROMA-output/denoised_func_data_nonaggr.nii.gz'])
    volBrain_file = ''.join([EPIpath,'rT1_brain_mask_FC.nii.gz'])
    regressors = np.array([])
    # set filename postfix for output image
    nR = 'aroma'

elif nuisanceReg == "HMParam":
    print("Applying Head Motion Param regressors")
    resting_file = ''.join([EPIpath,'4_epi.nii.gz'])
    volBrain_file = ''.join([EPIpath,'rT1_brain_mask_FC.nii.gz'])

    if numReg == 24:
        print(" -- 24 Head motion regressors")
        fname=''.join([EPIpath,'/HMPreg/motion12_regressors.npz'])
        m12reg = np.load(fname)
        print(sorted(m12reg.files))
        fname=''.join([EPIpath,'/HMPreg/motion_sq_regressors.npz'])
        m_sq_reg = np.load(fname)  
        print(sorted(m_sq_reg.files))
        #regressors = np.array([])
    elif numReg == 12:
        print(" -- 12 Head motion regressors")
        fname=''.join([EPIpath,'/HMPreg/motion12_regressors.npz'])
        m12reg = np.load(fname)

    else:
        print("Number of head motion regressors must be 12 or 24")

    # set filename postfix for output image
    nR = 'aroma'


END
}


##############################################################################

## PHYSIOLOGICAL REGRESSORS
echo "# =========================================================="
echo "# 5.3 APPLY REGRESSORS "
echo "# =========================================================="

if ${flags_EPI_NuisanceReg}; then
    if ${flags_NuisanceReg_AROMA}; then
        log "nuisanceReg AROMA"
        nuisanceReg="AROMA"
        configs_EPI_numReg=0
    elif ${flags_NuisanceReg_HeadParam}; then
        log "nuisanceReg HMParam"
        nuisanceReg="HMParam"        
    fi
else 
    nuisanceReg="none"
    configs_EPI_numReg=0
fi 



if ${flags_EPI_PhysiolReg}; then
    if ${flags_PhysiolReg_aCompCorr}; then  
        log "PhysiolReg - aCompCorr"
        physReg="aCompCorr"
        PhReg_path="${HMPpath}/aCompCorr" 
    elif ${flags_PhysiolReg_WM_CSF}; then
        log "PhysiolReg - Mean CSF & WM signal"
        physReg="Mn_WM_CSF"
        PhReg_path="${HMPpath}/PhysReg"
        configs_EPI_numPC=0
    fi 
else
    physReg="none"
    configs_EPI_numPC=0    
fi

if [[ ! ${flags_EPI_GS} ]]; then
    configs_EPI_numGS=0
fi

log "${EPIpath} ${nuisanceReg} ${configs_EPI_numPC} ${configs_EPI_numReg} ${configs_EPI_numGS} ${physReg}"
log "calling python script"
apply_reg ${EPIpath} ${nuisanceReg} ${configs_EPI_numPC} ${configs_EPI_numReg} ${configs_EPI_numGS} ${physReg}
