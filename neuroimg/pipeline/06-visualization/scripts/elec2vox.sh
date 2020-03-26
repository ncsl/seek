#!/bin/bash

while IFS=$'\t' read -r -a myArray; do
    x=$(printf "%.01f" "${myArray[3]}")
    y=$(printf "%.01f" "${myArray[4]}")
    z=$(printf "%.01f" "${myArray[5]}")
    z_=$(echo "256 - $z" | bc -l)
    seg=$(mri_info $SUBJECTS_DIR/$1/mri/aparc+aseg.mgz --voxel $x $y $z_)
    echo $seg
    segment=$(printf "%1.0f" "${seg}")
    while IFS=" " read -a line; do
        labelIndex=$(printf "%1.0f" ${line[0]})
        if [ "$segment" -eq "$labelIndex" ]; then
            echo $x $y $z_ "    " "${myArray[2]}" "    " ${line[1]} >> $SUBJECTS_DIR/$1/electrodes/electrodeLocations.txt
        fi
    done < LUT.txt
done < $SUBJECTS_DIR/$1/electrodes/electrodes.txt

