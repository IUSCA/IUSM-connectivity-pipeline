#!/bin/bash


# IU modules load
module unload python/2.7.16; module load python/3.6.8 
module load fsl/6.0.1; 
module load mricrogl
module load afni/18.3.03
module load ants
# module load singularity

################################################################################
################################################################################
## GLOBALS & dependencies

# where this package of scripts are
export EXEDIR=$(dirname "$(readlink -f "$0")")
# source bash funcs
# source activate py37_connPipeline
# source deactivate
source ${EXEDIR}/src/func/bash_funcs.sh
source ${EXEDIR}/config.sh

export python3_7="${EXEDIR}/../miniconda3/bin/python3.7"


# FSL
# set FSL env vars for fsl_sub.IU or fsl_sub.orig
if [[ -z ${FSLDIR} ]] ; then
	echoerr "FSLDIR not set"
	exit 1
fi


#################################################################################
# TODO -- List and loop through all subjects
export SUBJ="10692_1_AAK"
#export SUBJ="10692_1"
export T1path="${path2data}/${SUBJ}/${configs_T1}"
export DWIpath="${path2data}/${SUBJ}/${configs_DWI}"
export EPItemp="${path2data}/${SUBJ}/EPI"
#################################################################################


## main
main() {

start=`date +%s`

    # # get the whole call
    # cmdLineCall=$(echo "$0 $@")

    # # read in args
    # while (( $# > 1 )) ; do
    #     case "$1" in
    # 		-d | --dwi) shift
    # 			DWIinputnative="${1}" ; checkisfile ${DWIinputnative}
    # 			shift
    # 			;;
    # 		-o | --out) shift
    # 			OUTbasedir="${1}"
    # 			shift
    # 			;;
    #         -*)
    #             echo "ERROR: Unknown option '$1'"
    #             exit 1
    #             break
    #             ;;
    #         *)
    #             break
    #             ;;
    #      esac
    # done

    # shift "$((OPTIND-1))" # Shift off the options and optional

# ################################################################################
# ################################################################################
# ## check the argugments

    # # basic files needed
    # if [[ -z ${DWIinputnative} || -z ${OUTbasedir} ]] ; then
    # 	echoerr "missing required args: dwi img, bvecs, bvals, t1, out dir"
    # 	exit 1
    # fi

    # # get full path for out directroy
    # OUTbasedir=$(readlink -f ${OUTbasedir})

    # # cd into basedir
    # mkdir -vp ${OUTbasedir} || \
    # 	{ echoerr "could not make output dir. exitng" ; exit 1 ; }

    # echo ${cmdLineCall} > ${OUTbasedir}/cmdlinecall.txt

# ################################################################################
# ################################################################################
# ## run it

######################################################################################
log "# ############################ T1_PREPARE_A #####################################"

    if $T1_PREPARE_A; then

        cmd="${EXEDIR}/src/scripts/t1_prepare_A.sh" # -d ${PWD}/inputdata/dwi.nii.gz \
        echo $cmd
        eval $cmd
        exitcode=$?

        if [[ ${exitcode} -ne 0 ]] ; then
            echoerr "problem at T1_PREPARE_A. exiting."
            exit 1
        fi
    else 
        log "SKIP T1_PREPARE_A for subject $SUBJ"
    fi 

######################################################################################
log "# ############################ T1_PREPARE_B #####################################"


    if $T1_PREPARE_B; then

        if [[ -d "$T1path" ]]; then 

            cmd="${EXEDIR}/src/scripts/t1_prepare_B.sh" # -np ${numParcs} -d ${PWD}/inputdata/dwi.nii.gz \
            echo $cmd
            eval $cmd
            exitcode=$?

            if [[ ${exitcode} -ne 0 ]] ; then
                echoerr "problem at T1_PREPARE_B. exiting."
                exit 1
            fi
                    
        else
            echo "T1 directory doesn't exist; skipping subject $SUBJ"
        fi
    else
        log "SKIP T1_PREPARE_B for subject $SUBJ"
    fi 

######################################################################################
log "# ############################ fMRI_A ###########################################"


    if $fMRI_A; then

        if [[ -d "$T1path" ]]; then 

            cmd="${EXEDIR}/src/scripts/fMRI_A.sh"
            echo $cmd
            eval $cmd
            exitcode=$?

            if [[ ${exitcode} -ne 0 ]] ; then
                echoerr "problem at fMRI_A. exiting."
                exit 1
            fi
                    
        else
            echo "T1 directory doesn't exist; skipping subject $SUBJ"
        fi
    else
        log "SKIP fMRI_A for subject $SUBJ"
    fi 

######################################################################################
log "# ############################ fMRI_B ##########################################"

## Generates all the figures. Can still be called from Matlab for now...?


######################################################################################
log "# ############################ DWI_A ############################################"


    if $DWI_A; then

        if [[ -d "${DWIpath}" ]]; then 

            cmd="${EXEDIR}/src/scripts/DWI_A.sh"
            echo $cmd
            eval $cmd
            exitcode=$?

            if [[ ${exitcode} -ne 0 ]] ; then
                echoerr "problem at DWI_A. exiting."
                exit 1
            fi
                    
        else
            echo "T1 directory doesn't exist; skipping subject $SUBJ"
        fi
    else
        log "SKIP DWI_A for subject $SUBJ"
    fi 

# ################################################################################
# ################################################################################

    ## time it
    end=`date +%s`
    runtime=$((end-start))
    log "runtime: $runtime"

} # main

# ################################################################################
# ################################################################################
# ## run it

main "$@"


