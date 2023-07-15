#!/bin/bash

# ******************************************** USER INPUT ******************************************** #
jobname_suffix='neut2.5m_0.75z0_vtk_'

highfile='box_hr129921.nc'
highboxes=('HighT1_inflow0deg' 'HighT2_inflow0deg')

lowfile='box_lr129921.nc'
lowbox='Low'

offsetz=5  # see windtools/amr/postprocessing.py:to_vtk docstrings

# Optioanal inputs (requires a priori knowledge of how many boxes are there in total)
# To know amount of boxes: 
#     python ~/utilities/postprocess_amr_boxes2vtk.py -p . -f $highfile -g $highboxes
#     python ~/utilities/postprocess_amr_boxes2vtk.py -p . -f $lowfile  -g $lowbox
nNodes_low=20
nNodes_high=40  # Number of nodes per high box
itime_low=0
ftime_low=6001
itime_high=0
ftime_high=6001

# SLURM options
account='shellwind'
jobtime='4:00:00'

# **************************************************************************************************** #

# Deterine increments given number of nodes (already integers)
increment_low=$(((ftime_low-itime_low)/nNodes_low))
increment_high=$(((ftime_high-itime_high)/nNodes_high))   


# Make symlink if it doesn't exist
if ! [ -L postprocess_amr_boxes2vtk.py ]; then
    ln -s /home/rthedin/utilities/postprocess_amr_boxes2vtk.py
fi

# Get the current path. This script _needs_ to be launched from the case directory of interest
path=$(pwd -P)

echo -e "Current \$path is $path"

# Launch the high-res boxes on multiple nodes
for currgroup in ${highboxes[@]}; do
    curr_ftime=$itime_high
    for ((node=1;node<=nNodes_high;node++)); do
        # Get this node's starting and end time
        curr_itime=$curr_ftime
        curr_ftime=$((curr_itime+increment_high))
        # Special case for last chunk
        if [[ "$node" == "$nNodes_high" ]]; then
           curr_ftime=$ftime_high
        fi

        jobname=${jobname_suffix}${currgroup}_node${node}_ti${curr_itime}_tf${curr_ftime}
        echo -e "\nLaunching $jobname:\npostprocess_amr_boxes2vtk.py -p \$path -f $highfile -g $currgroup -offsetz $offsetz -itime $curr_itime -ftime $curr_ftime"
        sbatch -J $jobname -A $account -t $jobtime postprocess_amr_boxes2vtk.py -p $path -f $highfile -g $currgroup -offsetz $offsetz -itime $curr_itime -ftime $curr_ftime
    done
done


# Launch the low-res box on multiple nodes
#curr_ftime=$itime_low
#for ((node=1;node<=nNodes_low;node++)); do
#    # Get this node's starting and end time
#    curr_itime=$curr_ftime
#    curr_ftime=$((curr_itime+increment_low))
#    # Special case for last chunk
#    if [[ "$node" == "$nNodes_low" ]]; then
#        curr_ftime=$ftime_low
#    fi
#
#    jobname=${jobname_suffix}${lowbox}_node${node}_ti${curr_itime}_tf${curr_ftime}
#    echo -e "\nLaunching $jobname:\npostprocess_amr_boxes2vtk.py -p \$path -f $lowfile -g $lowbox -offsetz $offsetz -itime $curr_itime -ftime $curr_ftime"
#    sbatch -J $jobname -A $account -t $jobtime postprocess_amr_boxes2vtk.py -p $path -f $lowfile -g $lowbox -offsetz $offsetz -itime $curr_itime -ftime $curr_ftime
#done


