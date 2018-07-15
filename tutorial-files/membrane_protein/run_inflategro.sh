#!/bin/bash

# This script assumes you have already inflated the system with the InflateGRO command
# provided in the tutorial, and that you have further updated topol.top correctly to
# reflect the number of DPPC lipids that were deleted.
#
# This script assumes that the gmx binary is in the $PATH and that all necessary input
# files (coordinates, topology, .mdp) are in the current working directory.
#
# Confirm that you have achieved an appropriate area per lipid by inspecting the area_shrink*.dat
# files along the shrinking iterations.

# prevent accumulation of backup files like mdout.mdp
export GMX_MAXBACKUP=-1

# Run energy minimization on the inflated system
gmx grompp -f minim_inflategro.mdp -c system_inflated.gro -p topol.top -r system_inflated.gro -o system_inflated_em.tpr
gmx mdrun -deffnm system_inflated_em

# make molecules whole
echo 0 | gmx trjconv -s system_inflated_em.tpr -f system_inflated_em.gro -o tmp.gro -pbc mol
mv tmp.gro system_inflated_em.gro

# loop over 26 shrinking iterations
for curr in {1..26} 
do
    echo "########################################"
    echo "#"
    echo "# RUNNING SHRINKING ITERATION ${curr}..."
    echo "#"
    echo "########################################"

    prev=$((curr - 1))

    if [ $curr -eq 1 ]; then
        if [ ! -e system_inflated_em.gro ]; then
            echo "system_inflated_em.gro does not exist! Exiting."
            exit;
        fi
        # special file name if doing the first iteration
        perl inflategro.pl system_inflated_em.gro 0.95 DPPC 0 system_shrink${curr}.gro 5 area_shrink${curr}.dat
    else
        if [ ! -e system_shrink${prev}_em.gro ]; then
            echo "system_shrink${prev}_em.gro does not exist! Exiting."
            exit;
        fi
        # otherwise use minimized coordinates from previous iteration
        perl inflategro.pl system_shrink${prev}_em.gro 0.95 DPPC 0 system_shrink${curr}.gro 5 area_shrink${curr}.dat
    fi

    # run grompp and mdrun to carry out energy minimization
    gmx grompp -f minim_inflategro.mdp -c system_shrink${curr}.gro -r system_shrink${curr}.gro -p topol.top -o system_shrink${curr}_em.tpr
    gmx mdrun -deffnm system_shrink${curr}_em

    # make molecules whole
    gmx trjconv -s system_shrink${curr}_em.tpr -f system_shrink${curr}_em.gro -o tmp.gro -pbc mol
    mv tmp.gro system_shrink${curr}_em.gro

done

exit;
