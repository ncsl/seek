#!/bin/bash

# #PREP RECON
mkdir $SUBJECTS_DIR/$1/mri
# mkdir $SUBJECTS_DIR/$1/CT
mkdir $SUBJECTS_DIR/$1/electrodes
# mkdir $SUBJECTS_DIR/$1/Meshes
# mkdir $SUBJECTS_DIR/$1/Meshes/subcortical
mkdir $SUBJECTS_DIR/$1/mri/orig
# mkdir $SUBJECTS_DIR/$1/label/gyri
mkdir $SUBJECTS_DIR/$1/obj
mkdir $SUBJECTS_DIR/$1/rois
# # mv $SUBJECTS_DIR/$1/CT.nii $SUBJECTS_DIR/$1/CT/CT.nii
mv $SUBJECTS_DIR/$1/T1.nii $SUBJECTS_DIR/$1/mri/orig/T1.nii
mri_convert $SUBJECTS_DIR/$1/mri/orig/T1.nii $SUBJECTS_DIR/$1/mri/orig/001.mgz

./scripts/reconall.sh $1