
               
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
echo "# 3. Apply transformation to tissue and parcellation images."
echo "# =========================================================="

if [[ ! -e "${EPIpath}/T1_2_epi_dof6_bbr.mat" ]]; then  
    log "WARNING File ${EPIpath}/T1_2_epi_dof6_bbr.mat does not exist. Skipping further analysis"
    exit 1 
fi

#-------------------------------------------------------------------------#

# brain mask
fileIn="${T1path}/T1_brain_mask_filled.nii.gz"
fileRef="${EPIpath}/2_epi_meanvol_mask.nii.gz"
fileOut="${EPIpath}/rT1_brain_mask"
fileInit="${EPIpath}/T1_2_epi_dof6_bbr.mat"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd 

# WM mask
fileIn="${T1path}/T1_WM_mask.nii.gz"
fileOut="${EPIpath}/rT1_WM_mask"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd 

# Eroded WM mask
fileIn="${T1path}/T1_WM_mask_eroded.nii.gz"
fileOut="${EPIpath}/rT1_WM_mask_eroded"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd 

# CSF mask
fileIn="${T1path}/T1_CSF_mask.nii.gz"
fileOut="${EPIpath}/rT1_CSF_mask"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd 

# Eroded CSF mask
fileIn="${T1path}/T1_CSF_mask_eroded.nii.gz"
fileOut="${EPIpath}/rT1_CSF_mask_eroded"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd    

# CSF ventricle mask
fileIn="${T1path}/T1_CSFvent_mask_eroded.nii.gz"
fileOut="${EPIpath}/rT1_CSFvent_mask_eroded"
cmd="flirt -in ${fileIn} \
    -ref ${fileRef} \
    -out ${fileOut} \
    -applyxfm -init ${fileInit} \
    -interp nearestneighbour -nosearch"
log $cmd
eval $cmd  

#-------------------------------------------------------------------------#

# Probabilistic GM reg to fMRI spac
fileIn="${T1path}/T1_GM_mask.nii.gz"
fileOut="${EPIpath}/rT1_GM_mask_prob"
cmd="flirt -applyxfm \
-init ${fileInit} \
-interp spline \
-in ${fileIn} \
-ref ${fileRef} \
-out ${fileOut} -nosearch"
log $cmd
eval $cmd   

# binarize GM probability map
fileIn="${EPIpath}/rT1_GM_mask_prob.nii.gz"
fileOut="${EPIpath}/rT1_GM_mask"
cmd="fslmaths ${fileIn} -thr ${configs_EPI_GMprobthr} -bin ${fileOut}"
log $cmd
eval $cmd 

#-------------------------------------------------------------------------#
# Apllying T1 to EPI transformations to parcellations
for ((p=1; p<=numParcs; p++)); do  # exclude PARC0 - CSF - here

    parc="PARC$p"
    parc="${!parc}"
    pcort="PARC${p}pcort"
    pcort="${!pcort}"  
    pnodal="PARC${p}pnodal"  
    pnodal="${!pnodal}"                        

    echo "p is ${p} -- ${parc} parcellation -- pcort is -- ${pcort} -- pnodal is -- ${pnodal}"

    if [ ${pnodal} -eq 1 ]; then 
        
        # transformation from T1 to epi space
        fileIn="${T1path}/T1_GM_parc_${parc}_dil.nii.gz"
        fileOut="${EPIpath}/rT1_parc_${parc}.nii.gz"

        cmd="flirt -applyxfm -init ${fileInit} \
        -interp nearestneighbour \
        -in  ${fileIn} \
        -ref ${fileRef} \
        -out ${fileOut} -nosearch"
        log $cmd
        eval $cmd 

        # masking Shen with GM
        fileIn="${EPIpath}/rT1_parc_${parc}.nii.gz"                        
        fileOut="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"
        fileMul="${EPIpath}/rT1_GM_mask.nii.gz"

        cmd="fslmaths ${fileIn} \
        -mas ${fileMul} ${fileOut}"
        log $cmd
        eval $cmd                         
        
##================================================================================================
            ## UNFINISHED ##
            
        # removal of small clusters within ROIs
        fileIn="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"                        
        fileOut="${EPIpath}/rT1_GM_parc_${parc}_clean.nii.gz"   

        # #cmd="${EXEDIR}/src/scripts/get_largerst_clusters.sh ${fileIn} ${fileIn} ${configs_EPI_minVoxelsClust}"                     
        # cmd="python ${EXEDIR}/src/scripts/test_python_scripts.py"
        # log $cmd
        # eval $cmd 
##================================================================================================

    elif [ ${pnodal} -eq 0 ]; then # treat as an organizational parcellation to group nodes
        # Added by MDZ; 10/06/2015
        # transformation from T1 to epi space
        fileIn="${T1path}/T1_GM_parc_${parc}.nii.gz"
        fileOut="${EPIpath}/rT1_parc_${parc}.nii.gz"

        cmd="flirt -applyxfm -init ${fileInit} \
        -interp nearestneighbour \
        -in  ${fileIn} \
        -ref ${fileRef} \
        -out ${fileOut} -nosearch"
        log $cmd
        eval $cmd 

        # dilate parcellation
        fileIn="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"                        
        fileOut="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"

        cmd="fslmaths ${fileIn} \
        -dilD ${fileOut}"
        log $cmd
        eval $cmd        

        # masking parcellation with GM
        fileIn="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"                        
        fileOut="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"
        fileMul="${EPIpath}/rT1_GM_mask.nii.gz"

        cmd="fslmaths ${fileIn} \
        -mas ${fileMul} ${fileOut}"
        log $cmd
        eval $cmd      

##================================================================================================
            ## UNFINISHED ##
            
        # removal of small clusters within ROIs
        fileIn="${EPIpath}/rT1_GM_parc_${parc}.nii.gz"                        
        fileOut="${EPIpath}/rT1_GM_parc_${parc}_clean.nii.gz"   

        # #cmd="${EXEDIR}/src/scripts/get_largerst_clusters.sh ${fileIn} ${fileIn} ${configs_EPI_minVoxelsClust}"                     
        # cmd="python ${EXEDIR}/src/scripts/test_python_scripts.py"
        # log $cmd
        # eval $cmd 
##================================================================================================                                                               
    else
        log "WARNING the pnodal property is not specified for ${parc} parcellation.\
        Transformation to EPI not done"
    fi 

done 