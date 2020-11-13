"""
========================================
04b. Label Contacts Using SEEK Algorithm
========================================

This pipeline depends on the following functions:

    * python SEEK algorithm

"""

import os
import os.path as op
import sys
from pathlib import Path

from mne_bids import BIDSPath

sys.path.append("../../../")
from seek.pipeline.utils.fileutils import (BidsRoot, BIDS_ROOT,
                                           _get_seek_config, _get_ieeg_bids_dir,
                                           _get_anat_bids_dir, _get_bids_basename, )

configfile: _get_seek_config()

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=config.get('center_id'))
subject_wildcard = "{subject}"
SESSION = 'presurgery'

# initialize directories that we access
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_CT_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "CT"
FSOUT_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "elecs"
FSOUT_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "acpc"
FSOUT_SURF_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "surf"
FSOUT_GYRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / 'label' / 'gyri'

FIGDIR = os.path.join(FSOUT_ELECS_FOLDER, "figs")

# directory paths to output Raw BIDS related files
BIDS_PRESURG_IEEG_DIR = _get_ieeg_bids_dir(bids_root.bids_root, subject_wildcard, session=SESSION)
BIDS_PRESURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root, subject_wildcard, session='presurgery')
premri_fs_bids_fname = _get_bids_basename(subject_wildcard, session='presurgery',
                                          space='fs', imgtype='T1w', ext='nii')

t1_fs_fpath = os.path.join(BIDS_PRESURG_ANAT_DIR, premri_fs_bids_fname)

# the output electrode/coordsystem file(s)
#################### Centroid coordinates ####################
# manually labeled
manual_coordsystem_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='manual',
    suffix='coordsystem', extension='.json').basename
manual_electrodes_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='manual',
    suffix='electrodes', extension='.tsv').basename
manual_electrodes_json = manual_electrodes_fname.replace('.tsv', '.json')

# output of SEEK estimation algorithm
seek_coordsystem_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='seek',
    suffix='coordsystem', extension='.json').basename
seek_electrodes_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='seek',
    suffix='electrodes', extension='.tsv').basename
seek_electrodes_json = seek_electrodes_fname.replace('.tsv', '.json')

#################### Sub workflows ####################
subworkflow reconstruction_workflow:
    workdir:
           "../02-reconstruction/"
    snakefile:
             "../02-reconstruction/reconstruction.smk"
    configfile:
              _get_seek_config()

subworkflow contact_localization_workflow:
    workdir:
           "../04-contact_localization/"
    snakefile:
             "../04-contact_localization/contact_localization.smk"
    configfile:
              _get_seek_config()

rule label_electrode_anatomy:
    input:
         # bids_coordsystem_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, seek_coordsystem_fname), subject=subjects),
         # bids_electrodes_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, seek_electrodes_fname), subject=subjects),
         # bids_electrodes_json_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, seek_electrodes_json), subject=subjects),
         manual_coordsystem_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, manual_coordsystem_fname), subject=subjects),
         manual_bids_electrodes_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, manual_electrodes_fname), subject=subjects),
         manual_bids_electrodes_json_file=expand(op.join(BIDS_PRESURG_IEEG_DIR, manual_electrodes_json),
                                                 subject=subjects),
    params:
          bids_root=bids_root.bids_root,
    output:
          report=report('figcontactsanat.png', caption='report/figcontacts.rst', category='Contact Localization')
    shell:
         "echo 'done';"
         "bids-validator {params.bids_root};"
         "touch figcontactsanat.png {output};"

"""
Convert an annotation file into multiple label files or into a segmentation 'volume'. 
It can also create a border overlay.

# This version of mri_annotation2label uses the coarse labels from the Desikan-Killiany Atlas, unless
    # atlas_surf is 'destrieux', in which case the more detailed labels are used
"""
rule convert_gyri_annotations_to_labels:
    input:
         recon_success_file=reconstruction_workflow(os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")),
    params:
          subject_id=subject_wildcard,
          hemisphere='lh',
          surf_atlas_suffix_destrieux='--a2009s',
          surf_atlas_suffix_dk='',
          FS_GYRI_DIR=FSOUT_GYRI_FOLDER,
    shell:
         "mri_annotation2label " \
         "--subject {params.subject_id} " \
         "--hemi lh " \
         "--surface pial {params.surf_atlas_suffix_destrieux} " \
         "--outdir {FS_GYRI_DIR};"
         "mri_annotation2label " \
         "--subject {params.subject_id} " \
         "--hemi lh " \
         "--surface pial {params.surf_atlas_suffix_dk} " \
         "--outdir {FS_GYRI_DIR};"
         "mri_annotation2label " \
         "--subject {params.subject_id} " \
         "--hemi rh " \
         "--surface pial {params.surf_atlas_suffix_destrieux} " \
         "--outdir {FS_GYRI_DIR};"
         "mri_annotation2label " \
         "--subject {params.subject_id} " \
         "--hemi rh " \
         "--surface pial {params.surf_atlas_suffix_dk} " \
         "--outdir {FS_GYRI_DIR};"

rule apply_anatomicalatlas_to_electrodes:
    input:
         fs_lut_fpath=os.path.join('./', '"FreeSurferColorLUT.txt"),
         clustered_center_points_file=contact_localization_workflow(op.join(FSOUT_ELECS_FOLDER, seek_electrodes_fname)),
         T1_NIFTI_IMG=reconstruction_workflow(t1_fs_fpath),
         mri_coordsystem_fpath=contact_localization_workflow(op.join(FSOUT_ELECS_FOLDER, seek_coordsystem_fname)),
         params:
         FSPATIENT_DIR=str(FSPATIENT_SUBJECT_DIR),
         output:
         bids_electrodes_tsv_file=op.join(BIDS_PRESURG_IEEG_DIR, seek_electrodes_fname),
         bids_electrodes_json_file=op.join(BIDS_PRESURG_IEEG_DIR, seek_electrodes_json),
         mri_coordsystem_fpath=op.join(BIDS_PRESURG_IEEG_DIR, seek_coordsystem_fname),
         shell:
         "echo 'Applying anatomical atlas to electrode points...';"
         "python ./apply_anat_to_electrodes.py " \
             "{input.clustered_center_points_file} " \
             "{output.bids_electrodes_tsv_file} " \
             "{params.FSPATIENT_DIR} " \
             "{input.fs_lut_fpath} " \
             "{input.T1_NIFTI_IMG};"
         "cp {input.mri_coordsystem_fpath} {output.mri_coordsystem_fpath};"

rule apply_anatomicalatlas_to_manual_electrodes:
    input:
         fs_lut_fpath=os.path.join('./', "FreeSurferColorLUT.txt"),
         mri_bids_electrodes_tsv_file=contact_localization_workflow(
             op.join(FSOUT_ELECS_FOLDER, manual_electrodes_fname)),
         T1_NIFTI_IMG=reconstruction_workflow(t1_fs_fpath),
         mri_coordsystem_fpath=contact_localization_workflow(op.join(FSOUT_ELECS_FOLDER, manual_coordsystem_fname)),
    params:
          FSPATIENT_DIR=str(FSPATIENT_SUBJECT_DIR),
    output:
          bids_electrodes_tsv_file=op.join(BIDS_PRESURG_IEEG_DIR, manual_electrodes_fname),
          bids_electrodes_json_file=op.join(BIDS_PRESURG_IEEG_DIR, manual_electrodes_json),
          mri_coordsystem_fpath=op.join(BIDS_PRESURG_IEEG_DIR, manual_coordsystem_fname),
    shell:
         "echo 'Applying anatomical atlas to electrode points...';"
         "python ./apply_anat_to_electrodes.py " \
         "{input.mri_bids_electrodes_tsv_file} " \
         "{output.bids_electrodes_tsv_file} " \
         "{params.FSPATIENT_DIR} " \
         "{input.fs_lut_fpath} " \
         "{input.T1_NIFTI_IMG};"
         "cp {input.mri_coordsystem_fpath} {output.mri_coordsystem_fpath};"

         # rule apply_whitematteratlas_to_electrodes:
         #     input:
         #         WM_IMG_FPATH =  ,
         #     params:
         #         FSPATIENT_DIR = FSPATIENT_SUBJECT_DIR,
         #     output:
         #         bids_electrodes_file = os.path.join(electrodes_fname),
         #     shell:
         #         "echo 'Applying anatomical atlas to electrode points...';"
         #         "python ./apply_anat_to_electrodes.py " \
         #         "{input.clustered_center_points_file} " \
         #         "{input.clustered_center_voxels_file} " \
         #         "{output.bids_electrodes_file} " \
         #         "{params.FSPATIENT_DIR} " \
         #         "{input.WM_IMG_FPATH};"
