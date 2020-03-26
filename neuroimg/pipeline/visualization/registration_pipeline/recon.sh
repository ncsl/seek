#!/bin/bash

mkdir $SUBJECTS_DIR/$1/mri
mkdir $SUBJECTS_DIR/$1/mri/orig
mri_convert $SUBJECTS_DIR/$1/T1.nii $SUBJECTS_DIR/$1/mri/orig/001.mgz


if test -f "$SUBJECTS_DIR/$1/T2.nii"; then
    recon-all -subject $1 -T2 $SUBJECTS_DIR/$1/T2.nii -T2pial -parallel -openmp $(nproc) -all
else
    recon-all -subject $1 -parallel -openmp $(nproc) -all
fi
mkdir $SUBJECTS_DIR/$1/electrodes
mkdir $SUBJECTS_DIR/$1/obj
mkdir $SUBJECTS_DIR/$1/rois

