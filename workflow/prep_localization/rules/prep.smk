"""
===============================================
01. Prep Reconstruction Workflow and BIDS Layout
===============================================

In this pipeline, we prep the reconstruction workflow
by putting MRI and CT data into the BIDS layout and
re-orient images to RAS with ACPC alignment.

We assume that there is only one set of dicoms for CT and MRI
data.

This pipeline depends on the following functions:

    * mrconvert
    * acpcdetect

from FreeSurfer6+, acpcdetect2.0. To create a DAG pipeline, run:

    snakemake --dag | dot -Tpdf > dag_pipeline_reconstruction.pdf

"""

# Authors: Adam Li <adam2392@gmail.com>
# License: GNU

import os
import sys
from pathlib import Path

# hack to run from this file folder
sys.path.append("../../../")
from seek.utils.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_config,
                                  _get_ct_bids_dir,
                                  _get_bids_basename)

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']

# get the actual file path to the config
configpath = Path(_get_seek_config()).parent

# get the freesurfer patient directory
subject_wildcard = "{subject}"
bids_root = BidsRoot(subject_wildcard,BIDS_ROOT(config['bids_root']),
    site_id=config['site_id'],subject_wildcard=subject_wildcard)

# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
RAW_CT_FOLDER = bids_root.get_rawct_dir()
FSOUT_CT_FOLDER = Path(bids_root.get_freesurfer_patient_dir()) / "CT"
BIDS_PRESURG_CT_DIR = _get_ct_bids_dir(bids_root.bids_root,subject_wildcard,session='presurgery')

# original native files
ct_native_bids_fname = _get_bids_basename(subject_wildcard,
    session='presurgery',space='orig',
    imgtype='CT',ext='nii')

# output files
# raw T1/CT output
ct_output = os.path.join(BIDS_PRESURG_CT_DIR,ct_native_bids_fname)

logger.info('In prep localization workflow.')

rule prep:
    input:
        CT_bids_fname=expand(ct_output,subject=subjects),
    params:
        bids_root=bids_root.bids_root,
    log:
        expand("logs/prep_localization.{subject}.log",subject=subjects),
    output:
        report=report('fig1.png',caption='report/figprep.rst',category='Prep')
    shell:
        "echo 'done';"
        "bids-validator {params.bids_root};"
        "touch fig1.png {output};"

"""
Rule for prepping fs_recon by converting dicoms -> NIFTI images.

For more information, see BIDS specification.
"""

rule convert_dicom_to_nifti_ct:
    params:
        CT_FOLDER=RAW_CT_FOLDER,
        bids_root=bids_root.bids_root,
    log:
        "logs/prep_localization.{subject}.log"
    container:
        freesurfer_dockerurl
    output:
        CT_bids_fname=os.path.join(BIDS_PRESURG_CT_DIR,ct_native_bids_fname),
    shell:
        "mrconvert {params.CT_FOLDER} {output.CT_bids_fname};"
