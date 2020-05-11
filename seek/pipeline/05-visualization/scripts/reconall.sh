#!/bin/bash
cd $SUBJECTS_DIR/$1
recon-all -subjid $1 -openmp $(nproc) -all

