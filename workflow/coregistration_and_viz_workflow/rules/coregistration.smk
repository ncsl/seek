"""
======================================================
03. Coregistration (FSL-Flirt) from CT to T1 MRI
======================================================

This pipeline depends on the following functions:

    * flirt
    * mri_convert

from FreeSurfer6+, FSL.
"""

import os
import sys
from pathlib import Path

from mne_bids import BIDSPath

sys.path.append("../../../")
from seek.pipeline.utils.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_config,
                                           _get_bids_basename, _get_ct_bids_dir,
                                           _get_anat_bids_dir)

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']
fsl_dockerurl = config['fsl_docker']

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=config.get('center_id')
                     )
subject_wildcard = "{subject}"

# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_CT_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "CT"

BIDS_PRESURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root, subject_wildcard, session='presurgery')
BIDS_PRESURG_CT_DIR = _get_ct_bids_dir(bids_root.bids_root, subject_wildcard, session='presurgery')
ct_bids_fname = _get_bids_basename(subject_wildcard, session='presurgery',
                                   imgtype='CT', space='orig', ext='nii')

# MRI that FS uses
premri_fs_bids_fname = _get_bids_basename(subject_wildcard, session='presurgery',
                                          space='fs', imgtype='T1w', ext='nii')
t1_fs_fpath = os.path.join(BIDS_PRESURG_ANAT_DIR, premri_fs_bids_fname)

# CT mapped to FS in BIDs folders
ctint1_fs_bids_fname = _get_bids_basename(subject_wildcard,
                                          session='presurgery',
                                          space='fs',
                                          imgtype='CT', ext='nii')

from_id = 'CT'  # post implant CT
to_id = 'fs'  # freesurfer's T1w
kind = 'xfm'
ct_to_t1wfs_transform_fname = BIDSPath(subject=subject_wildcard,
                                       session='presurgery',
                                       space='fs').basename + \
                              f"_from-{from_id}_to-{to_id}_mode-image_{kind}.mat"
ct_tot1_fs_output = os.path.join(BIDS_PRESURG_CT_DIR, ctint1_fs_bids_fname)
ct_tot1_fs_map = os.path.join(BIDS_PRESURG_CT_DIR, ct_to_t1wfs_transform_fname)

print('In coregistration workflow.')

# subworkflow prep_workflow:
#     workdir:
#            "../01-prep/"
#     snakefile:
#              "../01-prep/prep.smk"
#     configfile:
#               _get_seek_config()

subworkflow reconstruction_workflow:
    workdir:
          "../02-reconstruction/"
    snakefile:
           "./02-reconstruction/reconstruction.smk"
    configfile:
              _get_seek_config()

# First rule
rule coregister_ct_and_T1w_images:
    input:
         # FLIRT FSL OUTPUT COREGISTRATION
         CT_IN_T1_NIFTI_IMG_ORIG=expand(os.path.join(FSOUT_CT_FOLDER, ctint1_fs_bids_fname), subject=subjects),
         # mapping matrix for CT to T1
         MAPPING_FILE=expand(os.path.join(FSOUT_CT_FOLDER, ct_to_t1wfs_transform_fname),
                             subject=subjects),
         # MAPPED BRAIN MASK TO CT SPACE
         brainmask_inct_file=expand(os.path.join(FSOUT_CT_FOLDER, "brainmask_inct.nii.gz"),
                                    subject=subjects),
         ct_in_fs_img=expand(ct_tot1_fs_output, subject=subjects),
         ct_in_fs_map=expand(ct_tot1_fs_map, subject=subjects),
    output:
          report=report('figct.png', caption='report/figprep.rst', category='Coregistration')
    shell:
         "echo 'done';"
         "touch figct.png {output};"

rule prep_ct_for_coregistration:
    input:
         CT_NIFTI_IMG=os.path.join(BIDS_PRESURG_CT_DIR, ct_bids_fname),
    params:
          CTDIR=str(FSOUT_CT_FOLDER),
    container:
             freesurfer_dockerurl
    output:
          CT_NIFTI_IMG=os.path.join(FSOUT_CT_FOLDER, ct_bids_fname),
    shell:
         "mri_convert {input.CT_NIFTI_IMG} {output.CT_NIFTI_IMG};"

"""
Rule for coregistering .nifit images -> .nifti for T1 space using Flirt in FSL.

E.g. useful for CT, and DTI images to be coregistered
"""
rule coregister_ct_to_t1wfs:
    input:
         PREMRI_NIFTI_IMG_MGZ=t1_fs_fpath,
         CT_NIFTI_IMG_MGZ=os.path.join(FSOUT_CT_FOLDER, ct_bids_fname),
    container:
             fsl_dockerurl
    output:
          # mapped image from CT -> MRI
          CT_IN_PRE_NIFTI_IMG_ORIGgz=os.path.join(FSOUT_CT_FOLDER, ctint1_fs_bids_fname + ".gz"),
          # mapping matrix for post to pre in T1
          MAPPING_FILE_ORIG=os.path.join(FSOUT_CT_FOLDER, ct_to_t1wfs_transform_fname),
          ct_tot1_fs_map=ct_tot1_fs_map,
    shell:
         "flirt -in {input.CT_NIFTI_IMG_MGZ} \
                             -ref {input.PREMRI_NIFTI_IMG_MGZ} \
                             -omat {output.MAPPING_FILE_ORIG} \
                             -out {output.CT_IN_PRE_NIFTI_IMG_ORIGgz};"
         "cp {output.MAPPING_FILE_ORIG} {output.ct_tot1_fs_map};"

rule convert_ctgz_to_nifti:
    input:
         CT_IN_PRE_NIFTI_IMG_ORIGgz=os.path.join(FSOUT_CT_FOLDER, ctint1_fs_bids_fname + ".gz"),
    container:
             fsl_dockerurl,
    output:
          CT_IN_PRE_NIFTI_IMG=os.path.join(FSOUT_CT_FOLDER, ctint1_fs_bids_fname),
          ct_tot1_fs_output=ct_tot1_fs_output,
    shell:
         "mrconvert {output.CT_IN_PRE_NIFTI_IMG_ORIGgz} {output.CT_IN_PRE_NIFTI_IMG};"
         "cp {output.CT_IN_PRE_NIFTI_IMG} {output.ct_tot1_fs_output};"

"""
Rule to map the brain mask over to the CT space.
"""
rule map_brainmask_to_ct:
    input:
         brainmask_file=os.path.join(FSOUT_MRI_FOLDER, "brainmask.nii.gz"),
         CT_NIFTI_IMG=os.path.join(FSOUT_CT_FOLDER, ct_bids_fname),
         # mapping matrix for post to pre in T1
         MAPPING_FILE_ORIG=os.path.join(FSOUT_CT_FOLDER, ct_to_t1wfs_transform_fname),
    container:
             fsl_dockerurl
    output:
          # mapping matrix for post to pre in T1
          brainmask_inct_file=os.path.join(FSOUT_CT_FOLDER, "brainmask_inct.nii.gz"),
    shell:
         "flirt -in {input.brainmask_file} \
                             -ref {input.CT_NIFTI_IMG} \
                             -applyxfm -init {input.MAPPING_FILE_ORIG} \
                             -out {output.brainmask_inct_file};"
