#!/bin/bash

## To run this script type: source red_startup_script.sh
# Do no run with ./ it will not have the same effect.

# First load python (required for FSL)
module unload python
module load python/3.6.8

# Load latest FSL version 6.0.1
module unload fsl
module load fsl/6.0.1

# Load afni 
module unload afni
module load afni/18.3.03
# It will try to load an older version of libpng than that loaded by FSL, but I am not sure if this will cause any issues. (2019.08.08 ejc)

#load matlab
module unload matlab
module load matlab/2019a

#load ants
module unload ants
module load ants

#load mrtrix
module load mrtrix/3.0
