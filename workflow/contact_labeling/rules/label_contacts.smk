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
from seek.utils.fileutils import (BidsRoot, BIDS_ROOT, DEFAULT_SESSION,
                                  _get_seek_config, _get_ieeg_bids_dir,
                                  _get_anat_bids_dir, _get_bids_basename, )

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']
fsl_dockerurl = config['fsl_docker']
seek_dockerurl = config['seek_docker']

# get the freesurfer patient directory
subject_wildcard = "{subject}"
coordframe_wildcard = "{coord_frame}"
coord_frames = ['mni', 'tkras']

bids_root = BidsRoot(subject_wildcard,BIDS_ROOT(config['bids_root']),
    site_id=config['site_id'],subject_wildcard=subject_wildcard)

# initialize directories that we access
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir()
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_CT_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "CT"
FSOUT_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "elecs"
FSOUT_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "acpc"
FSOUT_SURF_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "surf"
FSOUT_GYRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / 'label' / 'gyri'

FIGDIR = os.path.join(FSOUT_ELECS_FOLDER,"figs")

# directory paths to output Raw BIDS related files
BIDS_PRESURG_IEEG_DIR = _get_ieeg_bids_dir(bids_root.bids_root,subject_wildcard,session=DEFAULT_SESSION)
BIDS_PRESURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root,subject_wildcard,session=DEFAULT_SESSION)
premri_fs_bids_fname = _get_bids_basename(subject_wildcard,session='presurgery',
    space='fs',imgtype='T1w',ext='nii')

t1_fs_fpath = os.path.join(BIDS_PRESURG_ANAT_DIR,premri_fs_bids_fname)

# atlas image volumes to apply anatomical labels to
atlas_img_fpaths = [
    op.join(FSOUT_MRI_FOLDER,'aparc+aseg.mgz'),
    op.join(FSOUT_MRI_FOLDER,'aparc.a2009s+aseg.mgz'),
    op.join(FSOUT_MRI_FOLDER,'wmparc.mgz'),
]

# the output electrode/coordsystem file(s)
#################### Centroid coordinates ####################
# manually labeled
manual_coordsystem_fname = BIDSPath(
    subject=subject_wildcard,session=DEFAULT_SESSION,
    # processing='manual',
    space='mri',
    suffix='coordsystem',extension='.json').basename
manual_electrodes_fname = BIDSPath(
    subject=subject_wildcard,session=DEFAULT_SESSION,
    # processing='manual',
    space='mri',
    suffix='electrodes',extension='.tsv').basename
manual_electrodes_json = manual_electrodes_fname.replace('.tsv','.json')

# mni, tkras
other_coordsystem_fname = BIDSPath(
    subject=subject_wildcard,session=DEFAULT_SESSION,
    # processing='manual',
    space=coordframe_wildcard,
    suffix='coordsystem',extension='.json').basename
other_electrodes_fname = BIDSPath(
    subject=subject_wildcard,session=DEFAULT_SESSION,
    # processing='manual',
    space=coordframe_wildcard,
    suffix='electrodes',extension='.tsv').basename

rule label_electrode_anatomy:
    input:
        manual_coordsystem_file=expand(op.join(BIDS_PRESURG_IEEG_DIR,manual_coordsystem_fname),
            subject=subjects,coord_frame=coord_frames),
        manual_bids_electrodes_file=expand(op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_fname),
            subject=subjects,coord_frame=coord_frames),
        manual_bids_electrodes_json_file=expand(op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_json),
            subject=subjects,coord_frame=coord_frames),
        other_coordsystem_fpath=expand(op.join(BIDS_PRESURG_IEEG_DIR,other_coordsystem_fname),
            subject=subjects,coord_frame=coord_frames),
        other_electrodes_tsv_file=expand(op.join(BIDS_PRESURG_IEEG_DIR,other_electrodes_fname),
            subject=subjects,coord_frame=coord_frames),
    params:
        bids_root=bids_root.bids_root,
    log: expand("logs/label-contacts.{subject}.log",subject=subjects)
    container:
        seek_dockerurl
    output:
        report=report('figcontactsanat.png',caption='report/figcontacts.rst',category='Contact Localization')
    shell:
        "echo 'done';"
        "bids-validator {params.bids_root};"
        "touch figcontactsanat.png {output};"


rule convert_gyri_annotations_to_labels:
    """
    Convert an annotation file into multiple label files or into a segmentation 'volume'. 
    It can also create a border overlay.

    # This version of mri_annotation2label uses the coarse labels from the Desikan-Killiany Atlas, unless
        # atlas_surf is 'destrieux', in which case the more detailed labels are used
    """
    input:
        recon_success_file=os.path.join(FSPATIENT_SUBJECT_DIR,"{subject}_recon_success.txt"),
    params:
        subject_id=subject_wildcard,
        hemisphere='lh',
        surf_atlas_suffix_destrieux='--a2009s',
        surf_atlas_suffix_dk='',
        FS_GYRI_DIR=FSOUT_GYRI_FOLDER,
    log: "logs/label-contacts.{subject}.log",
    container:
        freesurfer_dockerurl
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


rule label_desikanatlas_anatomy_on_electrodes:
    """Labels anatomy for electrodes.tsv files.
    
    Assumes electrodes are localized according to the
    FreeSurfer ``T1.mgz`` file (i.e. FreeSurfer T1 space).
    
    """
    input:
        mri_bids_electrodes_tsv_file=(
            op.join(FSOUT_ELECS_FOLDER,manual_electrodes_fname)),
        mri_coordsystem_fpath=op.join(FSOUT_ELECS_FOLDER,manual_coordsystem_fname),
        atlas_img=expand('{image_fpath}',image_fpath=atlas_img_fpaths),
    params:
        bids_root=str(bids_root.bids_root),
    log: "logs/label-contacts.{subject}.log",
    container:
        seek_dockerurl
    output:
        bids_electrodes_tsv_file=op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_fname),
        bids_electrodes_json_file=op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_json),
        bids_mri_coordsystem_fpath=op.join(BIDS_PRESURG_IEEG_DIR,manual_coordsystem_fname),
    shell:
        "echo 'Applying anatomical atlas to electrode points...';"
        "cp {input.mri_coordsystem_fpath} {output.bids_mri_coordsystem_fpath};"
        "python ./labeling_scripts/label_anatomy.py " \
        "   {input.mri_bids_electrodes_tsv_file} " \
        "   {input.atlas_img}"
        "   {output.bids_electrodes_tsv_file} " \
        "   {params.bids_root};"
        "echo 'Copying coordsystem.json file from derivatives folder to BIDS root...';"


rule save_mni_coordinates_electrodes:
    input:
        bids_electrodes_tsv_file=op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_fname),
        bids_electrodes_json_file=op.join(BIDS_PRESURG_IEEG_DIR,manual_electrodes_json),
        bids_mri_coordsystem_fpath=op.join(BIDS_PRESURG_IEEG_DIR,manual_coordsystem_fname),
    params:
        bids_root=str(bids_root.bids_root),
        coord_frame=coordframe_wildcard,
    log: "logs/label-contacts.{subject}.log",
    container:
        seek_dockerurl
    output:
        other_coordsystem_fpath=op.join(BIDS_PRESURG_IEEG_DIR,other_coordsystem_fname),
        other_electrodes_tsv_file=op.join(BIDS_PRESURG_IEEG_DIR,other_electrodes_fname),
    shell:
        "echo 'Converting coordinate systems of electrode points...';"
        "cp {input.mri_coordsystem_fpath} {output.bids_mri_coordsystem_fpath};"
        "python ./labeling_scripts/convert_coordsystem.py " \
        "   {input.mri_bids_electrodes_tsv_file} " \
        "   {input.atlas_img}"
        "   {output.bids_electrodes_tsv_file} " \
        "   {params.bids_root};"
