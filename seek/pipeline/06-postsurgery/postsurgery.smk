"""
===============================================
06. Postsurgery T1 MRI registration and workflow
================================================

This pipeline depends on the following functions:

    * mrconvert
    * flirt

from FreeSurfer6+, FSL.
"""

import os
import sys
from pathlib import Path

from mne_bids import make_bids_basename

sys.path.append("../../../")
from seek.pipeline.utils.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_config,
                                           _get_anat_bids_dir, _get_bids_basename)

configfile: _get_seek_config()

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=config.get('center_id'))
subject_wildcard = "{subject}"

SESSION = 'postsurgery'

# initialize directories that we access in this snakemake
RAW_CT_FOLDER = bids_root.get_rawct_dir(subject_wildcard)
RAW_POSTMRI_FOLDER = bids_root.get_postmri_dir(subject_wildcard)
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_POSTMRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "postsurgerymri"
FSOUT_POSTMRI_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "postsurgerymri"

BIDS_POSTSURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root, subject_wildcard, session=SESSION)
BIDS_PRESURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root, subject_wildcard, session='presurgery')

###########################
# ORIGINAL IMAGES
# original space
postmri_native_bids_fname = _get_bids_basename(subject_wildcard,
                                               session=SESSION, space='orig',
                                               imgtype='T1w', ext='nii')
premri_native_bids_fname = _get_bids_basename(subject_wildcard,
                                              session='presurgery', space='orig',
                                              imgtype='T1w', ext='nii')

# robust fov file
premri_robustfov_native_bids_fname = _get_bids_basename(subject_wildcard,
                                                        session='presurgery',
                                                        space='orig',
                                                        processing='robustfov',
                                                        imgtype='T1w', ext='nii')
postmri_robustfov_native_bids_fname = _get_bids_basename(subject_wildcard,
                                                         session=SESSION,
                                                         space='orig',
                                                         processing='robustfov',
                                                         imgtype='T1w', ext='nii')

# after ACPC
premri_bids_fname = _get_bids_basename(subject_wildcard,
                                       session='presurgery',
                                       space='ACPC',
                                       imgtype='T1w', ext='nii')
postmri_bids_fname = _get_bids_basename(subject_wildcard,
                                        session=SESSION,
                                        space='ACPC',
                                        imgtype='T1w', ext='nii')

# after FreeSurfer
premri_fs_bids_fname = _get_bids_basename(subject_wildcard,
                                        session='presurgery',
                                        space='fs',
                                        imgtype='T1w', ext='nii')

###########################
# COREGISTERED IMAGES
# Coregistration filenames
postinpre_bids_fname = _get_bids_basename(subject_wildcard,
                                          session=SESSION,
                                          space='T1w',
                                          imgtype='T1w', ext='nii')

from_id = 'postT1w'  # post implant CT
to_id = 'preT1w'  # freesurfer's T1w
kind = 'xfm'
pre_to_post_transform_fname = make_bids_basename(subject=subject_wildcard,
                                                 session=SESSION,
                                                 space='T1w') + \
                              f"_from-{from_id}_to-{to_id}_mode-image_{kind}.mat"

# after ACPC
postinpre_acpc_bids_fname = _get_bids_basename(subject_wildcard,
                                               session=SESSION,
                                               space='T1wACPC',
                                               imgtype='T1w', ext='nii')

from_id = 'postT1w'  # post implant CT
to_id = 'preT1w'  # freesurfer's T1w
kind = 'xfm'
pre_to_post_acpc_transform_fname = make_bids_basename(subject=subject_wildcard,
                                                      session=SESSION,
                                                      space='T1wACPC') + \
                                   f"_from-{from_id}_to-{to_id}_mode-image_{kind}.mat"

# after FreeSurfer
postinpre_fs_bids_fname = _get_bids_basename(subject_wildcard,
                                             session=SESSION,
                                             space='fs',
                                             imgtype='T1w', ext='nii')

from_id = 'postT1w'  # post implant CT
to_id = 'fs'  # freesurfer's T1w
kind = 'xfm'
pre_to_post_fs_transform_fname = make_bids_basename(subject=subject_wildcard,
                                                    session=SESSION,
                                                    space='fs') + \
                                 f"_from-{from_id}_to-{to_id}_mode-image_{kind}.mat"

# output files
# raw original T1w output
t1_output = os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname)
t1_acpc_output = os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_bids_fname)

# raw post to pre T1 image and map
post_tot1_output = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                                postinpre_bids_fname)
post_tot1_map = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                             pre_to_post_transform_fname)

# post to pre T1acpc image and map
post_tot1_acpc_output = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                                     postinpre_acpc_bids_fname)
post_tot1_acpc_map = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                                  pre_to_post_acpc_transform_fname)

# post to pre T1 FS image and map
post_tot1_fs_output = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                                   postinpre_fs_bids_fname)
post_tot1_fs_map = os.path.join(BIDS_POSTSURG_ANAT_DIR,
                                pre_to_post_fs_transform_fname)

subworkflow prep_workflow:
    workdir:
           "../01-prep/"
    snakefile:
             "../01-prep/prep.smk"
    configfile:
              _get_seek_config()

# First rule
rule postprep:
    input:
         MRI_NIFTI_IMG=expand(t1_output, subject=subjects),
         post_acpc_output=expand(t1_acpc_output, subject=subjects),
         post_tot1_output=expand(post_tot1_output, subject=subjects),
         post_tot1_map=expand(post_tot1_map, subject=subjects),
         post_tot1_acpc_output=expand(post_tot1_acpc_output, subject=subjects),
         post_tot1_acpc_map=expand(post_tot1_acpc_map, subject=subjects),
         # post_tot1_fs_output=expand(post_tot1_fs_output, subject=subjects),
         # post_tot1_fs_map=expand(post_tot1_fs_map, subject=subjects),
    params:
          bids_root=bids_root.bids_root,
    output:
          report=report('figpost.png', caption='report/figpost.rst', category='Prep')
    shell:
         "echo 'done';"
         "bids-validator {params.bids_root};"
         "touch fig1.png {output};"

"""
Rule for prepping fs_recon by converting dicoms -> NIFTI images.

They are reoriented into 'LAS' orientation. For more information, see
BIDS specification.
"""
rule convert_dicom_to_bids:
    params:
          MRI_FOLDER=RAW_POSTMRI_FOLDER,
          bids_root=bids_root.bids_root,
    output:
          MRI_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname),
    shell:
         "mrconvert {params.MRI_FOLDER} {output.MRI_bids_fname};"

"""
Apply robust FOV.
"""
rule postt1w_compute_robust_fov:
    input:
         MRI_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname),
    output:
          MRI_bids_fname_gz=os.path.join(FSOUT_POSTMRI_ACPC_FOLDER, postmri_robustfov_native_bids_fname) + '.gz',
          MRI_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_robustfov_native_bids_fname),
    container:
             "docker://neuroseek/seek",
    shell:
         "echo 'robustfov -i {input.MRI_bids_fname} -r {output.MRI_bids_fname_gz}';"
         'robustfov -i {input.MRI_bids_fname} -r {output.MRI_bids_fname_gz};'  # -m roi2full.mat
         'mrconvert {output.MRI_bids_fname_gz} {output.MRI_bids_fname} --datatype uint16;'

"""
Rule for automatic ACPC alignment using acpcdetect software. 

Please check the output images to quality assure that the ACPC was properly
aligned.  
"""
rule postt1w_automatic_acpc_alignment:
    input:
         MRI_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_robustfov_native_bids_fname),
    params:
          anat_dir=str(BIDS_POSTSURG_ANAT_DIR),
          acpc_fs_dir=str(FSOUT_POSTMRI_ACPC_FOLDER),
    output:
          MRI_bids_fname_fscopy=os.path.join(FSOUT_POSTMRI_ACPC_FOLDER, postmri_bids_fname),
          MRI_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_bids_fname),
    shell:
         # create BIDS session directory and copy file there
         "echo 'acpcdetect -i {input.MRI_bids_fname} -center-AC -output-orient RAS;'"
         "echo {output.MRI_bids_fname};"
         "mkdir -p {params.acpc_fs_dir};"
         "cp {input.MRI_bids_fname} {output.MRI_bids_fname_fscopy};"
         # run acpc auto detection
         "acpcdetect -i {output.MRI_bids_fname_fscopy} -center-AC -output-orient RAS;"
         "cp {output.MRI_bids_fname_fscopy} {output.MRI_bids_fname};"

"""
Rule for coregistering .nifit images -> .nifti for T1 space using Flirt in FSL.

E.g. useful for CT, and DTI images to be coregistered
"""
rule coregister_t1w_post_to_pre:
    input:
         pre_bids_fname=os.path.join(BIDS_PRESURG_ANAT_DIR, premri_native_bids_fname),
         post_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname),
    params:
          FSOUT_POSTMRI_FOLDER=str(FSOUT_POSTMRI_FOLDER),
    output:
          # mapped image from CT -> MRI
          POST_IN_PRE_NIFTI_IMG_ORIGgz=os.path.join(FSOUT_POSTMRI_FOLDER, postinpre_bids_fname + ".gz"),
          POST_IN_PRE_NIFTI_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, postinpre_bids_fname),
          # mapping matrix for post to pre in T1
          MAPPING_FILE_ORIG=os.path.join(FSOUT_POSTMRI_FOLDER, pre_to_post_transform_fname),
          MAPPING_FILE_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, pre_to_post_transform_fname),
    shell:
         "mkdir --parents {params.FSOUT_POSTMRI_FOLDER};"
         "flirt -in {input.post_bids_fname} \
                             -ref {input.pre_bids_fname} \
                             -omat {output.MAPPING_FILE_ORIG} \
                             -out {output.POST_IN_PRE_NIFTI_IMG_ORIGgz};"
         "mrconvert {output.POST_IN_PRE_NIFTI_IMG_ORIGgz} {output.POST_IN_PRE_NIFTI_BIDS};"
         "cp {output.MAPPING_FILE_ORIG} {output.MAPPING_FILE_BIDS};"

"""
Rule for coregistering .nifit images -> .nifti for T1 space using Flirt in FSL.

E.g. useful for CT, and DTI images to be coregistered
"""
rule coregister_t1w_post_to_preacpc:
    input:
         pre_bids_fname=os.path.join(BIDS_PRESURG_ANAT_DIR, premri_bids_fname),
         post_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname),
    params:
          FSOUT_POSTMRI_FOLDER=str(FSOUT_POSTMRI_FOLDER),
    output:
          # mapped image from CT -> MRI
          POST_IN_PRE_NIFTI_IMG_ORIGgz=os.path.join(FSOUT_POSTMRI_FOLDER, postinpre_acpc_bids_fname + ".gz"),
          POST_IN_PRE_NIFTI_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, postinpre_acpc_bids_fname),
          # mapping matrix for post to pre in T1
          MAPPING_FILE_ORIG=os.path.join(FSOUT_POSTMRI_FOLDER, pre_to_post_acpc_transform_fname),
          MAPPING_FILE_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, pre_to_post_acpc_transform_fname),
    shell:
         "mkdir --parents {params.FSOUT_POSTMRI_FOLDER};"
         "flirt -in {input.post_bids_fname} \
                             -ref {input.pre_bids_fname} \
                             -omat {output.MAPPING_FILE_ORIG} \
                             -out {output.POST_IN_PRE_NIFTI_IMG_ORIGgz};"
         "mrconvert {output.POST_IN_PRE_NIFTI_IMG_ORIGgz} {output.POST_IN_PRE_NIFTI_BIDS};"
         "cp {output.MAPPING_FILE_ORIG} {output.MAPPING_FILE_BIDS};"

"""
Rule for coregistering .nifit images -> .nifti for T1 space using Flirt in FSL.

E.g. useful for CT, and DTI images to be coregistered
"""
rule coregister_t1w_post_to_FSt1w:
    input:
         pre_bids_fname=os.path.join(BIDS_PRESURG_ANAT_DIR, premri_fs_bids_fname),
         post_bids_fname=os.path.join(BIDS_POSTSURG_ANAT_DIR, postmri_native_bids_fname),
    params:
          FSOUT_POSTMRI_FOLDER=str(FSOUT_POSTMRI_FOLDER),
    output:
          # mapped image from CT -> MRI
          POST_IN_PRE_NIFTI_IMG_ORIGgz=os.path.join(FSOUT_POSTMRI_FOLDER, postinpre_fs_bids_fname + ".gz"),
          POST_IN_PRE_NIFTI_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, postinpre_fs_bids_fname),
          # mapping matrix for post to pre in T1
          MAPPING_FILE_ORIG=os.path.join(FSOUT_POSTMRI_FOLDER, pre_to_post_fs_transform_fname),
          MAPPING_FILE_BIDS=os.path.join(BIDS_POSTSURG_ANAT_DIR, pre_to_post_fs_transform_fname),
    shell:
         "mkdir --parents {params.FSOUT_POSTMRI_FOLDER};"
         "flirt -in {input.post_bids_fname} \
                             -ref {input.pre_bids_fname} \
                             -omat {output.MAPPING_FILE_ORIG} \
                             -out {output.POST_IN_PRE_NIFTI_IMG_ORIGgz};"
         "mrconvert {output.POST_IN_PRE_NIFTI_IMG_ORIGgz} {output.POST_IN_PRE_NIFTI_BIDS};"
         "cp {output.MAPPING_FILE_ORIG} {output.MAPPING_FILE_BIDS};"
