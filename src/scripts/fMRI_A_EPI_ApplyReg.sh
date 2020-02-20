
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
EPIpath="$1" nuisanceReg="$2" config_param="$3" numReg="$4" numGS="$5" physReg="$6" ${python3_7} - <<END
import os
import numpy as np
import nibabel as nib
from scipy import stats

print("inside Python script")


EPIpath=os.environ['EPIpath']
nuisanceReg=os.environ['nuisanceReg']
print("nuisanceReg",nuisanceReg)
config_param=int(os.environ['config_param'])
print("config_param",config_param)
numReg=int(os.environ['numReg'])
print("numReg",numReg)
numGS=int(os.environ['numGS'])
print("numGS",numGS)
physReg=os.environ['physReg']
print("physReg",physReg)

PhReg_path = ''.join([EPIpath,'/HMPreg/'])
PhReg_path = ''.join([PhReg_path,physReg])
print("PhReg_path ",PhReg_path )


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
        regressors = np.concatenate((m12reg['motion'],m12reg['motion_deriv']),axis=1)
        regressors = np.concatenate((regressors,m_sq_reg['motion_sq']),axis=1)
        regressors = np.concatenate((regressors,m_sq_reg['motion_deriv_sq']),axis=1)
        print("regressors shape ",regressors.shape)
    elif numReg == 12:
        print(" -- 12 Head motion regressors")
        fname=''.join([EPIpath,'/HMPreg/motion12_regressors.npz'])
        m12reg = np.load(fname)
        print(sorted(m12reg.files))
        regressors = np.concatenate((m12reg['motion'],m12reg['motion_deriv']),axis=1)
        print("regressors shape ",regressors.shape)

    else:
        print("Number of head motion regressors must be 12 or 24")

if numGS > 0:
    fname = ''.join([PhReg_path,'/dataGS.npz'])
    dataGS = np.load(fname) 
    if numGS == 1:
        regressors = np.concatenate((regressors,dataGS['GSavg']),axis=1)
        print("   -- 1 global signal regressor ")
        print("regressors shape ",regressors.shape)
    if numGS == 2:
        regressors = np.concatenate((regressors,dataGS['GSavg']),axis=1)
        regressors = np.concatenate((regressors,dataGS['GSderiv']),axis=1)
        print("   -- 2 global signal regressor ")
        print("regressors shape ",regressors.shape)
    if numGS == 4:
        regressors = np.concatenate((regressors,dataGS['GSavg']),axis=1)
        regressors = np.concatenate((regressors,dataGS['GSavg_sq']),axis=1)
        regressors = np.concatenate((regressors,dataGS['GSderiv']),axis=1)
        regressors = np.concatenate((regressors,dataGS['GSderiv_sq']),axis=1)
        print("   -- 4 global signal regressor ")   
        print("regressors shape ",regressors.shape)             


if physReg == "aCompCorr":
    fname = ''.join([PhReg_path,'/dataPCA_WM-CSF.npz'])
    numphys = np.load(fname) 
    print("-- aCompCor PC of WM & CSF regressors")
    zRegressMat = [];
    if config_param > 5:
        CSFpca = numphys['CSFpca'].T
        WMpca = numphys['WMpca'].T
        print(" -- Applying all levels of PCA removal")
        for ic in range(6):
            if ic == 0:
                zRegressMat.append(stats.zscore(regressors,axis=0));
            else:
                regMat = np.concatenate((regressors,CSFpca[:,:ic]),axis=1)
                regMat = np.concatenate((regMat,WMpca[:,:ic]),axis=1)
                zRegressMat.append(stats.zscore(regMat,axis=0));
                print("num elements in zRegressMat is ",len(zRegressMat))

    elif 0 < config_param < 6:
        print(" -- Writing prespecified removal of %d components" % config_param)
        zRegressMat.append(stats.zscore(regressors,axis=0));



elif physReg == "PhysReg":
    fname = ''.join([PhReg_path,'/dataMnRg_WM-CSF.npz'])
    numphys = np.load(fname) 
    if config_param == 2:
        regressors = np.concatenate((regressors,numphys['CSFavg']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMavg']),axis=1)
        print("   -- 2 physiological regressors ")
        print("regressors shape ",regressors.shape)
    if config_param == 4:
        regressors = np.concatenate((regressors,numphys['CSFavg']),axis=1)
        regressors = np.concatenate((regressors,numphys['CSFderiv']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMavg']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMderiv']),axis=1)
        print("   -- 4 physiological regressors")
        print("regressors shape ",regressors.shape)
    if config_param == 8:
        regressors = np.concatenate((regressors,numphys['CSFavg']),axis=1)
        regressors = np.concatenate((regressors,numphys['CSFavg_sq']),axis=1)
        regressors = np.concatenate((regressors,numphys['CSFderiv']),axis=1)
        regressors = np.concatenate((regressors,numphys['CSFderiv_sq']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMavg']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMavg_sq']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMderiv']),axis=1)
        regressors = np.concatenate((regressors,numphys['WMderiv_sq']),axis=1)
        print("   -- 8 physiological regressors ")   
        print("regressors shape ",regressors.shape) 
    zRegressMat = [];
    zRegressMat.append(stats.zscore(regressors,axis=0));

    # print(len(zRegressMat))
    # print(type(zRegressMat))
    # print(type(zRegressMat[0]))


END
}


##############################################################################

## PHYSIOLOGICAL REGRESSORS
echo "# =========================================================="
echo "# 5.3 APPLY REGRESSORS "
echo "# =========================================================="

# if ${flags_EPI_NuisanceReg}; then
    if ${flags_NuisanceReg_AROMA}; then
        log "nuisanceReg AROMA"
        nuisanceReg="AROMA"
        configs_EPI_numReg=0
    elif ${flags_NuisanceReg_HeadParam}; then
        log "nuisanceReg HMParam"
        nuisanceReg="HMParam"        
    fi
# else 
#     nuisanceReg="none"
#     configs_EPI_numReg=0
# fi 



if ${flags_EPI_PhysiolReg}; then
    if ${flags_PhysiolReg_aCompCorr}; then  
        log "PhysiolReg - aCompCorr"
        physReg="aCompCorr"
        #configs_EPI_numPhys=0
        config_param=${configs_EPI_numPC}
    elif ${flags_PhysiolReg_WM_CSF}; then
        log "PhysiolReg - Mean CSF & WM signal"
        physReg="PhysReg" #"Mn_WM_CSF"
        config_param=${configs_EPI_numPhys}
        #configs_EPI_numPC=0
    fi 
else
    physReg="none"
    #configs_EPI_numPC=0   
    #configs_EPI_numPhys=0 
    config_param=0
fi

if [[ ! ${flags_EPI_GS} ]]; then
    configs_EPI_numGS=0
fi

log "${EPIpath} ${nuisanceReg} ${config_param} ${configs_EPI_numReg} ${configs_EPI_numGS} ${physReg}"
log "calling python script"
apply_reg ${EPIpath} ${nuisanceReg} ${config_param} ${configs_EPI_numReg} ${configs_EPI_numGS} ${physReg}
