"""
=============================================
04. Contact Localization Using SEEK Algorithm
=============================================

This pipeline depends on the following functions:

    * python SEEK algorithm

"""

import os
import os.path as op
import sys
from pathlib import Path

from mne_bids import BIDSPath, make_bids_folders

sys.path.append("../../../")
from seek.pipeline.utils.fileutils import (BidsRoot, BIDS_ROOT,
                                           _get_seek_config,
                                           _get_anat_bids_dir)

configfile: _get_seek_config()

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=config.get('center_id'))
subject_wildcard = "{subject}"
SESSION = 'presurgery'

# the pattern of the filename that gets output for each
# subject's Electrode Annotation Via FieldTrip Toolbox
elec_mat_fname_pattern = f'{subject}CH_elec_initialize.mat'


# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_CT_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "CT"
FSOUT_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "elecs"
FSOUT_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "acpc"
FSOUT_SURF_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "surf"
FSOUT_GYRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / 'label' / 'gyri'

FIGDIR = os.path.join(FSOUT_ELECS_FOLDER, "figs")


def _get_bids_basename(subject, acquisition):
    """Wildcard function to get bids_basename."""
    bids_fname = BIDSPath(subject, acquisition=acquisition,
                          suffix="electrodes", extension=".tsv")
    return bids_fname


def get_channels_tsv_fpath(bids_root, subject_id):
    channels_tsv_fpath = BIDSPath(root=bids_root, subject=subject_id,
                                  suffix="channels", extension=".tsv").fpath
    return channels_tsv_fpath

# input files
brainmask_ct_fpath = os.path.join(FSOUT_CT_FOLDER, "brainmask_inct.nii.gz")
ct_nifti_fpath = os.path.join(FSOUT_CT_FOLDER, "CT.nii")

# the path to the subject/session
data_path = make_bids_folders(subject=subject_wildcard, session=SESSION,
                              bids_root=str(bids_root.bids_root), make_dir=False,
                              overwrite=False, verbose=False)

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

#################### Voxel coordinates ####################
# manually labeled
voxel_manual_electrodes_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='manual',
    suffix='electrodes', extension='.tsv').basename

# output of SEEK estimation algorithm
voxel_electrodes_fname = BIDSPath(
    subject=subject_wildcard, session=SESSION, processing='seek',
    suffix='electrodes', extension='.tsv').basename

#################### Sub workflows ####################
subworkflow reconstruction_workflow:
    workdir:
           "../02-reconstruction/"
    snakefile:
             "../02-reconstruction/reconstruction.smk"
    configfile:
              _get_seek_config()

subworkflow coregistration_workflow:
    workdir:
           "../03-coregistration/"
    snakefile:
             "../03-coregistration/coregistration.smk"
    configfile:
              _get_seek_config()

#################### Rules for this pipeline ####################
# First rule
rule localize_electrode_coordinates:
    input:
         # figure_file=expand(os.path.join(FIGDIR, "summary_pca_elecs.png"),
         #                    subject=subjects),
         # elecs3d_figure_fpath=expand(os.path.join(FIGDIR, "electrodes_in_voxel_space.png"),
         #                             subject=subjects),
         bids_coordsystem_file=expand(op.join(FSOUT_ELECS_FOLDER, seek_coordsystem_fname), subject=subjects),
         bids_electrodes_file=expand(op.join(FSOUT_ELECS_FOLDER, seek_electrodes_fname), subject=subjects),
         bids_electrodes_json_file=expand(op.join(FSOUT_ELECS_FOLDER, seek_electrodes_json), subject=subjects),
         manual_coordsystem_file=expand(op.join(FSOUT_ELECS_FOLDER, manual_coordsystem_fname), subject=subjects),
         manual_bids_electrodes_file=expand(op.join(FSOUT_ELECS_FOLDER, manual_electrodes_fname), subject=subjects),
         manual_bids_electrodes_json_file=expand(op.join(FSOUT_ELECS_FOLDER, manual_electrodes_json), subject=subjects),
    params:
          bids_root=bids_root.bids_root,
    output:
          report=report('figcontacts.png', caption='report/figcontacts.rst', category='Contact Localization')
    shell:
         "echo 'done';"
         "bids-validator {params.bids_root};"
         "touch figcontacts.png {output};"

"""
Rule for conversion .mat -> .txt

Converts output of fieldtrip toolbox .mat files named accordingly into txt files. 
Note: This file should contain at least two contacts per electrode.
"""
rule convert_CT_eleccoords_to_txt:
    # input:
    #     clustered_center_points_mat = os.path.join(FSOUT_ELECS_FOLDER,
    #                                                # '{subject}_elec_f.mat'
    #                                                '{subject}CH_elec_initialize.mat'),
    params:
          elecs_dir=FSOUT_ELECS_FOLDER,
    output:
          clustered_center_points=os.path.join(FSOUT_ELECS_FOLDER, '{subject}_elecxyz_inct.txt'),
    shell:
         "python ./convert_to_txt.py {params.elecs_dir} {output.clustered_center_points};"


"""Semi-automated algorithm to find electrode contacts on CT image.

Requires at 2 end-point contacts on each electrode.
"""
rule find_electrodes_on_CT:
    input:
         CT_NIFTI_IMG=reconstruction_workflow(ct_nifti_fpath),
         brainmask_inct_file=coregistration_workflow(brainmask_ct_fpath),
         # list of channel points (at least 2 per electrode)
         electrode_initialization_file=os.path.join(FSOUT_ELECS_FOLDER, elec_mat_fname_pattern),
    params:
          FSPATIENT_DIR=FSPATIENT_SUBJECT_DIR.as_posix(),
    output:
          clustered_center_points=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(ct_electrodes_fname)),
          clustered_center_voxels=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(ct_voxel_electrodes_fname)),
          binarized_ct_volume=os.path.join(FSOUT_CT_FOLDER, "{subject}_binarized_ct.nii.gz"),
          elecs3d_figure_fpath=os.path.join(FIGDIR, "electrodes_in_voxel_space.png"),
          ct_coordsystem_fname=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(ct_coordsystem_fname)),
    shell:
         "echo 'RUNNING CLUSTERING ALGORITHM';"
         "python ./localize.py " \
         "{input.CT_NIFTI_IMG} " \
         "{input.brainmask_inct_file} " \
         "{input.electrode_initialization_file} " \
         "{output.clustered_center_points} " \
         "{output.clustered_center_voxels} " \
         "{output.binarized_ct_volume} " \
         "{params.FSPATIENT_DIR} " \
         "{output.elecs3d_figure_fpath};"


"""
Rule for image space conversion CT -> T1

applying flirt rigid registration affine transformation to the xyz coordinates of the localized
contacts in CT space. This will convert them into the space of the T1 image.
"""
rule convert_seeklabels_to_native_T1:
    input:
         CT_NIFTI_IMG=os.path.join(FSOUT_CT_FOLDER, "CT.nii"),
         MRI_NIFTI_IMG=os.path.join(FSOUT_MRI_FOLDER, "T1.nii"),
         # DEPENDENCY ON RECONSTRUCTION WORKFLOW
         # mapping matrix for post to pre in T1
         MAPPING_FILE=coregistration_workflow(os.path.join(FSOUT_CT_FOLDER, "fsl_ct-to-t1_omat.txt")),
         ct_xyzcoords_fpath=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(ct_electrodes_fname)),
    output:
          mri_xyzcoords_fpath=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(seek_electrodes_fname)),
          mri_coordsystem_fpath=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(seek_coordsystem_fname)),
    shell:
         "python ./convert_coordspace.py " \
         "{input.CT_NIFTI_IMG} " \
         "{input.MRI_NIFTI_IMG} " \
         "{input.MAPPING_FILE} " \
         "{input.ct_xyzcoords_fpath} " \
         "{output.mri_xyzcoords_fpath};"

rule convert_manuallabels_to_native_T1:
    input:
         CT_NIFTI_IMG=os.path.join(FSOUT_CT_FOLDER, "CT.nii"),
         MRI_NIFTI_IMG=os.path.join(FSOUT_MRI_FOLDER, "T1.nii"),
         # DEPENDENCY ON RECONSTRUCTION WORKFLOW
         # mapping matrix for post to pre in T1
         MAPPING_FILE=coregistration_workflow(os.path.join(FSOUT_CT_FOLDER, "fsl_ct-to-t1_omat.txt")),
         manual_ct_coords_fpath=os.path.join(FSOUT_ELECS_FOLDER, '{subject}_elecxyz_inct.txt'),
    output:
          mri_electrode_coords_file=os.path.join(FSOUT_ELECS_FOLDER, manual_electrodes_fname),
          mri_coordsystem_fpath=os.path.join(FSOUT_ELECS_FOLDER, manual_coordsystem_fname),
    shell:
         "python ./convert_coordspace.py " \
         "{input.CT_NIFTI_IMG} " \
         "{input.MRI_NIFTI_IMG} " \
         "{input.MAPPING_FILE} " \
         "{input.manual_ct_coords_fpath} " \
         "{output.mri_electrode_coords_file};"

"""
Rule to plot:
- parcellated nodes in space with the surface shown transparently
- surface shown transparent with different regions colored
- parcellated nodes in space with the contacts.xyz (centers) plotted
"""
rule visualize_results:
    input:
         ct_xyzcoords_fpath=os.path.join(FSOUT_ELECS_FOLDER, os.path.basename(ct_electrodes_fname)),
         # list of channel points (at least 2 per electrode)
         manual_electrode_file=os.path.join(FSOUT_ELECS_FOLDER, '{subject}_elecxyz_inct.txt'),
         # DEPENDENCY ON RECONSTRUCTION WORKFLOW
         ctfile=reconstruction_workflow(os.path.join(FSOUT_CT_FOLDER, "CT.nii")),
    params:
          fsdir=FSPATIENT_SUBJECT_DIR,
    output:
          pcafigure_file=os.path.join(FIGDIR, "summary_pca_elecs.png"),
          l2figure_file=os.path.join(FIGDIR, "summary_euclidean_distance_errors.png"),
    shell:
         "echo 'RUNNING VISUALIZATION CHECKS';"
         "python ./visualize_results.py {input.ct_xyzcoords_fpath} " \
         "{input.manual_electrode_file} " \
         "{input.ctfile} " \
         "{params.fsdir} " \
         "{output.pcafigure_file} {output.l2figure_file};"

"""Copy over the channels.tsv file to FreeSurfer."""
rule copy_channnelstsv_to_freesurfer:
    input:
         channels_tsv_fpath=lambda wildcards: get_channels_tsv_fpath(str(bids_root.bids_root),
                                                                     subject_id=wildcards.subject),
    output:
          channels_tsv_fpath=os.path.join(FSOUT_ELECS_FOLDER,
                                          os.path.basename(get_channels_tsv_fpath(str(bids_root.bids_root),
                                                                                  subject_id=subject_wildcard))),
    shell:
         "cp {input.channels_tsv_fpath} {output.channels_tsv_fpath};"
