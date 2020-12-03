"""
============================================
02. Reconstruction Workflow using FreeSurfer
============================================

In this pipeline, we prep the reconstruction workflow
by putting MRI and CT data into the BIDS layout and
re-orient images to RAS with ACPC alignment.

We assume that there is only one set of dicoms for CT and MRI
data.

This pipeline depends on the following functions:

    * mrconvert
    * acpcdetect

from FreeSurfer6+, acpcdetect2.0.
"""

import os
import sys
from pathlib import Path

sys.path.append("../../../")
from seek.pipeline.utils.fileutils import SCRIPTS_UTIL_DIR
from seek.pipeline.utils.fileutils import (BidsRoot, BIDS_ROOT,
                                           _get_seek_config,
                                           _get_anat_bids_dir, _get_bids_basename,
                                           _get_subject_center)

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']
fsl_dockerurl = config['fsl_docker']

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=_get_subject_center(subjects, centers, subject))
subject_wildcard = "{subject}"

# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_DIR = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FSOUT_MRI_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "mri"
FSOUT_CT_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "CT"
FSOUT_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "elecs"
FSOUT_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "acpc"
FSOUT_SURF_FOLDER = Path(FSPATIENT_SUBJECT_DIR) / "surf"

BIDS_PRESURG_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root, subject_wildcard, session='presurgery')
premri_bids_fname = _get_bids_basename(subject_wildcard, session='presurgery',
                                       space='ACPC', imgtype='T1w', ext='nii')

print('In reconstruction workflow.')

subworkflow prep_workflow:
    workdir:
           "../"
    snakefile:
             "./rules/prep.smk"
    configfile:
              _get_seek_config()

# output files of this pipeline
rhpial_asc_fpath = os.path.join(FSOUT_SURF_FOLDER, "ascii", "rh.pial.asc")
lhpial_asc_fpath = os.path.join(FSOUT_SURF_FOLDER, "ascii", "lh.pial.asc")
brainmask_nii_fpath = os.path.join(FSOUT_MRI_FOLDER, "brainmask.nii.gz")

premri_fs_bids_fname = _get_bids_basename(subject_wildcard, session='presurgery',
                                          space='fs', imgtype='T1w', ext='nii')
t1_fs_fpath = os.path.join(BIDS_PRESURG_ANAT_DIR, premri_fs_bids_fname)

# First rule
rule recon_all_output:
    input:
         outsuccess_file=expand(os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt"),
                                subject=subjects),
         subcort_success_flag_file=expand(os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_subcort_success.txt"),
                                          subject=subjects),
         brainmask_nifti=expand(brainmask_nii_fpath, subject=subjects),
         lhpial=expand(lhpial_asc_fpath, subject=subjects),
         rhpial=expand(rhpial_asc_fpath, subject=subjects),
         t1_fs=expand(t1_fs_fpath, subject=subjects)
    output:
          report=report('figreconstruct.png', caption='report/figpost.rst', category='FreeSurfer')
    shell:
         "echo 'done';"
         "touch figreconstruct.png {output};"

"""
Rule for prepping fs_recon.

The purpose is to setup a directory file structure that plays nice w/ Freesurfer. If you are using
some other module it is recommended to just add directories in there.

Rule for converting .dicom -> .niftiPUT_DIR,
                            "lh.native.aparc.annot"),
        rhlabel=os.path.join(NATIVESPACE_OUTPUT_DIR,
                            "rh.native.aparc.annot"),

"""
rule prep_recon:
    input:
         MRI_NIFTI_IMG=prep_workflow(os.path.join(BIDS_PRESURG_ANAT_DIR, premri_bids_fname)),
    params:
          elecsdir=str(FSOUT_ELECS_FOLDER),
          acpcdir=str(FSOUT_ACPC_FOLDER),
          INPUT_FS_ORIG_DIR=os.path.join(FSOUT_MRI_FOLDER, "orig"),
    container:
        freesurfer_dockerurl
    output:
          MRI_MGZ_IMG=os.path.join(FSOUT_MRI_FOLDER, "orig", "001.mgz"),
    shell:
         "mkdir -p {params.elecsdir};"
         "mkdir -p {params.acpcdir};"
         "mkdir -p {params.INPUT_FS_ORIG_DIR};"
         "mri_convert {input.MRI_NIFTI_IMG} {output.MRI_MGZ_IMG};"

"""
Rule for reconstructions .nifti -> output files.

Since Freesurfer creates the directory on its own + snakemake does too,
I instead specify an output as a "temporary" flagger file that will let snakemake
know that reconstruction was completed.
"""
rule reconstruction:
    input:
         MRI_MGZ_IMG=os.path.join(FSOUT_MRI_FOLDER, "orig", "001.mgz"),
    params:
          patient=subject_wildcard,
          SUBJECTS_DIR=FS_DIR,
    container:
        freesurfer_dockerurl
    output:
          outsuccess_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")
    shell:
         "export SUBJECTS_DIR={params.SUBJECTS_DIR};" \
         "SUBJECTS_DIR={params.SUBJECTS_DIR};"
         # "recon-all " \
         # "-cw256 " \
         #     # "-i {input.MRI_MGZ_IMG} " \
         # "-subject {params.patient} " \
         # "-all " \
         # "-parallel -openmp $(nproc); " \
         # "--bids-out;"
         "touch {output.outsuccess_file}"

"""
Rule for converting the pial surfaces to ascii data, so that it is readable by python/matlab.
"""
rule convert_pial_surface_files_ascii:
    input:
         recon_success_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")
    params:
          lhpial=os.path.join(FSOUT_SURF_FOLDER, "lh.pial"),
          rhpial=os.path.join(FSOUT_SURF_FOLDER, "rh.pial"),
    container:
        freesurfer_dockerurl
    output:
          lhpial=os.path.join(FSOUT_SURF_FOLDER, "ascii", "lh.pial.asc"),
          rhpial=os.path.join(FSOUT_SURF_FOLDER, "ascii", "rh.pial.asc"),
    shell:
         "mris_convert {params.lhpial} {output.lhpial};"
         "mris_convert {params.rhpial} {output.rhpial};"

"""
Rule for converting the pial surfaces to ascii data, so that it is readable by python/matlab.
"""
rule convert_FS_T1_to_nifti:
    input:
         recon_success_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")
    params:
          PREMRI_MGZ_IMG=os.path.join(FSOUT_MRI_FOLDER, "T1.mgz"),
    container:
        freesurfer_dockerurl
    output:
          PREMRI_NIFTI_IMG=os.path.join(FSOUT_MRI_FOLDER, "T1_fs_LIA.nii"),
          t1_fs_fpath=t1_fs_fpath
    shell:
         "mri_convert {params.PREMRI_MGZ_IMG} {output.PREMRI_NIFTI_IMG};"
         "mri_convert {params.PREMRI_MGZ_IMG} {output.t1_fs_fpath};"

"""
Rule for converting brainmask image volume from MGZ to Nifti.
"""
rule convert_brainmask_to_nifti:
    input:
         recon_success_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")
    params:
          brainmask_mgz=os.path.join(FSOUT_MRI_FOLDER, "brainmask.mgz")
    container:
        freesurfer_dockerurl
    output:
          brainmask_nifti=os.path.join(FSOUT_MRI_FOLDER, "brainmask.nii.gz")
    shell:
         "mri_convert {params.brainmask_mgz} {output.brainmask_nifti};"

"""
Rule for extracting the subcortical regions

- creates a new folder aseg2srf inside the fs_output data directory
- c/p this into our final result directory
"""
rule create_subcortical_volume:
    input:
         outsuccess_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt"),
         subcort_fs7_segmentation_success_flag_file=os.path.join(FSPATIENT_SUBJECT_DIR,
                                                                 f'{subject_wildcard}_fs7_subcort_segmentation_success.txt'),
    params:
          SUBJECTS_DIR=FS_DIR,
          patient=subject_wildcard,
          scripts_dir=SCRIPTS_UTIL_DIR,
    container:
         freesurfer_dockerurl
    output:
          subcort_success_flag_file=os.path.join(FSPATIENT_SUBJECT_DIR, f'{subject_wildcard}_subcort_success.txt'),
    shell:
         # generate subcortical region volume bounding surfaces
         "export SUBJECTS_DIR={params.SUBJECTS_DIR};"
         "SUBJECTS_DIR={params.SUBJECTS_DIR};"
         "{params.scripts_dir}/aseg2srf -s {params.patient};"
         "touch {output.subcort_success_flag_file};"

"""
Rule for FreeSurfer 7+, where one can segment the brainstem, hippocampus and thalamus further.

References:
- http://freesurfer.net/fswiki/ThalamicNuclei
- http://freesurfer.net/fswiki/HippocampalSubfieldsAndNucleiOfAmygdala
- http://freesurfer.net/fswiki/BrainstemSubstructures
"""
rule create_subcortical_segmentations:
    input:
         outsuccess_file=os.path.join(FSPATIENT_SUBJECT_DIR, "{subject}_recon_success.txt")
    params:
          SUBJECTS_DIR=FS_DIR,
          subject=subject_wildcard,
          scripts_dir=SCRIPTS_UTIL_DIR,
    container:
        freesurfer_dockerurl
    output:
#           output_files=os.path.join(FSPATIENT_SUBJECT_DIR, 'mri', 'ThalamicNuclei.v12.T1.volumes.txt'),
#             ThalamicNuclei.v12.T1.mgz: stores the discrete segmentation volumes at subvoxel resolution (0.5 mm).
# ThalamicNuclei.v12.T1.FSvoxelSpace.mgz:
          subcort_fs7_segmentation_success_flag_file=os.path.join(FSPATIENT_SUBJECT_DIR,
                                                                  f'{subject_wildcard}_fs7_subcort_segmentation_success.txt'),
    shell:
         "export SUBJECTS_DIR={params.SUBJECTS_DIR};"
         "SUBJECTS_DIR={params.SUBJECTS_DIR};"
         # Run these only after recon-all
         "segmentThalamicNuclei.sh {params.subject};"
         "segmentHA_T1.sh {params.subject};"
         "segmentBS.sh {params.subject};"
         "touch {output.subcort_fs7_segmentation_success_flag_file};"
