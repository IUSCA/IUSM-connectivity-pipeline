#!/bin/bash
#
#
#
###############################################################################
#
# ---------------------------- GET_READOUT ---------------------------------
# Obtain readout time from dicom data for distortion correction in FSL EDDY.
#
#                   Effective Echo Spacing (s) = 
#        1 / (BandwidthPerPixelPhaseEncode (Hz) * MatrixSizePhase)
#
#   BandwidthPerPixelPhaseEncode -> Siemens dicom tag (0019, 1028)
#   MatrixSizePhase -> size of image in phase encode direction; usually the
# first number in the field (0051, 100b), AcquisitionMatrixText

#               Total Readout Time (FSL definition)
# Time from the center of the first echo to the center of the last
# =(actual number of phase encoding lines - 1)* effective echo spacing
#
# Actual number of phase encoding lines = 
#                       Image matrix in phase direction / GRAPPA factor

# Original Matlab code - Evgeny Chumin, Indiana University School of Medicine, 2018
#                        John West, Indiana University School of Medicine, 2018
##
###############################################################################

# shopt -s nullglob # No-match globbing expands to null

# source ${EXEDIR}/src/func/bash_funcs.sh

###############################################################################

function largest_clusters() {
RedoutTime=$(python - "$1"<<END

END
)
}


###############################################################################

## Set paths and check for dicom direcotry
local fIn=$1
local fOut=$2
local threshold=$3

#echo "Calling Python script 
largest_clusters ${fIn}
echo "${RedoutTime}"


