"""
======================================================
07. Group Analysis Workflow using FreeSurfer to MNI152
======================================================

This pipeline depends on the following functions:

    * mri_cvs_register

from FreeSurfer6+
"""

# fsVox2RAS = np.array(
#             [[-1., 0., 0., 128.], [0., 0., 1., -128.], [0., -1., 0., 128.], [0., 0., 0., 1.]])

import os
import sys
from pathlib import Path

from mne_bids import make_bids_basename, make_bids_folders

sys.path.append("../../../")
from seek.pipeline.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_config,
                                     TEMPLATE_SUBJECT, _get_session_name)

configfile: _get_seek_config()

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']))
subject_wildcard = "{subject}"

# initialize directories that we access in this snakemake
RAW_CT_FOLDER = bids_root.get_rawct_dir(subject_wildcard)
RAW_MRI_FOLDER = bids_root.get_postmri_dir(subject_wildcard)
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_POSTMRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "postsurgerymri"
FSOUT_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "elecs"
FSOUT_CVS_DIR = Path(FSPATIENT_SUBJECT_DIR) / 'cvs'

session = _get_session_name(config)
acquisition = 'seeg'
# acquisition = _get_acquisition_name(config)

data_path = make_bids_folders(subject=subject_wildcard, session=session,
                              bids_root=str(bids_root.bids_root), make_dir=False,
                              overwrite=False, verbose=False)
manual_coordsystem_fname = make_bids_basename(
    subject=subject_wildcard, session=session, processing='manual',
    acquisition=acquisition, space='fs',
    suffix='coordsystem.json')
manual_electrodes_fname = make_bids_basename(
    subject=subject_wildcard, session=session, processing='manual',
    acquisition=acquisition, space='fs',
    suffix='electrodes.tsv')

cvs_electrodes_fname = make_bids_basename(
    subject=subject_wildcard, session=session, processing='manual',
    acquisition=acquisition, space='mni',
    suffix='electrodes.tsv', prefix=data_path)

TEMP_DIR = '/tmp'
TEMP_LABELS_TO_WARP_DIR = os.path.join(TEMP_DIR, 'labels_to_warp')
TEMP_WARPED_LABELS_DIR = os.path.join(TEMP_DIR, 'warped_labels')

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

# First rule
rule all:
    input:
         # morph file from {subject} to template subject
         morph_file=expand(os.path.join(FSOUT_CVS_DIR, f'combined_to{TEMPLATE_SUBJECT}_elreg_afteraseg-norm.tm3d'),
                           subject=config['patients']),
    params:
          bids_root=bids_root.bids_root,
    shell:
         "echo 'done';"
         "bids-validator {params.bids_root};"

"""
Rule for applying nonlinear warping subject to the MNI coordinate template 
using  FreeSurfer's mri_cvs_register (i.e. combined volume and surface) 
registration.
"""
rule cvs_register:
    input:
         PREMRI_NIFTI_IMG=reconstruction_workflow(os.path.join(FSOUT_MRI_FOLDER, "T1_fs_LIA.nii")),
    params:
          subject=subject_wildcard,
          template=TEMPLATE_SUBJECT,
          SUBJECTS_DIR=FS_DIR,
    output:
          morph_file=os.path.join(FSOUT_CVS_DIR, f'combined_to{TEMPLATE_SUBJECT}_elreg_afteraseg-norm.tm3d'),
    shell:
         "export SUBJECTS_DIR={params.SUBJECTS_DIR};"
         "echo 'Performing mri cvs nonlinear warping to MNI template...';"
         'mri_cvs_register --mov {params.subject} ' \
         '--template {params.template} ' \
         '--openmp 4 ' \
         '--nocleanup;'

         # rule prep_cvs_warp_of_electrodes:
         #     input:
         #         template_mri_mgz=os.path.join(FS_DIR, TEMPLATE_SUBJECT, 'mri', 'brain.mgz'),
         #         morph_file=os.path.join(FSOUT_CVS_DIR, f'combined_to{TEMPLATE_SUBJECT}_elreg_afteraseg-norm.tm3d'),
         #         fs_electrodes_fpath = contact_localization_workflow(os.path.join(FSOUT_ELECS_FOLDER, manual_electrodes_fname)),
         #     output:
         #         rh_label_fpath = os.path.join(TEMP_LABELS_TO_WARP_DIR, f'rh_{manual_electrodes_fname}'),
         #         lh_label_fpath = os.path.join(TEMP_LABELS_TO_WARP_DIR, f'lh_{manual_electrodes_fname}'),
         #         depth_label_fpath = os.path.join(TEMP_LABELS_TO_WARP_DIR, f'depth_{manual_electrodes_fname}'),
         #     shell:

         # rule apply_cvs_warp_to_depth_electrodes:
         #     input:
         #          template_mri_mgz=os.path.join(FS_DIR, TEMPLATE_SUBJECT, 'mri', 'brain.mgz'),
         #          morph_file=os.path.join(FSOUT_CVS_DIR, f'combined_to{TEMPLATE_SUBJECT}_elreg_afteraseg-norm.tm3d'),
         #          RAS_depth_label_fpath = os.path.join(TEMP_LABELS_TO_WARP_DIR, f'depth_{manual_electrodes_fname}'),
         #     output:
         #           nearest_point_coords_file=os.path.join(TEMP_WARPED_LABELS_DIR, cvs_electrodes_fname),
         #     shell:
         #          "echo 'Applying nonlinear warping to electrodes...;"
         #          'applyMorph --template {input.template_mri_mgz} ' \
         #          '--transform {input.morph_file} ' \
         #          'tract_point_list {input.RAS_depth_label_fpath} ' \
         #          '{output.nearest_point_coords_file} nearest;'

         # rule apply_cvs_warp_to_surface_electrodes:
         #     input:
         #          rh_label_fpath=os.path.join(TEMP_LABELS_TO_WARP_DIR, f'rh_{manual_electrodes_fname}'),
         #          lh_label_fpath=os.path.join(TEMP_LABELS_TO_WARP_DIR, f'lh_{manual_electrodes_fname}'),
         #     params:
         #           subject=subject_wildcard,
         #           template=TEMPLATE_SUBJECT,
         #           subject_dir=FSPATIENT_SUBJECT_DIR,
         #     output:
         #           rh_trglabel_fpath=os.path.join(TEMP_WARPED_LABELS_DIR, f'rh_{manual_electrodes_fname}'),
         #           lh_trglabel_fpath=os.path.join(TEMP_WARPED_LABELS_DIR, f'lh_{manual_electrodes_fname}'),
         #     shell:
         #          'mri_label2label --srclabel {input.rh_label_fpath} --srcsubject {params.subject} ' \
         #          '--trgsubject {params.template} ' \
         #          '--trglabel {output.rh_trglabel_fpath} ' \
         #          ' --regmethod surface --hemi rh ' \
         #          ' --trgsurf pial --paint 6 pial --sd {params.subject_dir};'
         #          'mri_label2label --srclabel {input.lh_label_fpath} --srcsubject {params.subject} ' \
         #          '--trgsubject {params.template} ' \
         #          '--trglabel {output.lh_trglabel_fpath} ' \
         #          ' --regmethod surface --hemi rh ' \
         #          ' --trgsurf pial --paint 6 pial --sd {params.subject_dir};'
