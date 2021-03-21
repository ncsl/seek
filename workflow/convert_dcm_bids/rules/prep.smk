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

from FreeSurfer6+, acpcdetect2.0.

Docker containers
-----------------
To run the workflow, there are a series of containers that we
depend on. They are containers for:

    - acpcdetectv2.0
    - FSLv6.0
    - FreeSurfer7.0 with MRtrix3v3.0.2

Notes
-----
To create a DAG pipeline, run:

    snakemake --dag | dot -Tpdf > dag_pipeline_reconstruction.pdf
"""

# Authors: Adam Li <adam2392@gmail.com>
# License: GNU

import os
import sys
from pathlib import Path

from mne_bids import BIDSPath

# hack to run from this file folder
sys.path.append("../../../")
from seek.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_config,
                            _get_anat_bids_dir, _get_ct_bids_dir,
                            _get_bids_basename)

configfile: _get_seek_config()

logger.debug('In prep workflow...')

# get the actual file path to the config
configpath = Path(_get_seek_config()).parent

# get the freesurfer patient directory
subject_wildcard = "{subject}"
session = config['session']
if session is None:
    session = 'presurgery'
bids_root = BidsRoot(subject_wildcard,BIDS_ROOT(config['bids_root']),
    site_id=config['site_id'],subject_wildcard=subject_wildcard)

freesurfer_dockerurl = config['freesurfer_docker']
acpcdetect_dockerurl = config['acpcdetect_docker']
fsl_dockerurl = config['fsl_docker']
seek_dockerurl = config['seek_docker']
heudiconv_dockerurl = config['heudiconv_docker']
validator_dockerurl = config['validator_docker']

# initialize directories that we access in this snakemake
sourcedata_path = bids_root.bids_root / 'sourcedata'
FS_DIR = bids_root.freesurfer_dir
RAW_CT_FOLDER = bids_root.get_rawct_dir()
RAW_MRI_FOLDER = bids_root.get_premri_dir()

BIDS_ANAT_DIR = _get_anat_bids_dir(bids_root.bids_root,subject_wildcard,session=session)
BIDS_CT_DIR = _get_ct_bids_dir(bids_root.bids_root,subject_wildcard,session=session)

# original native files
ct_native_bids_fname = BIDSPath(
    subject=subject_wildcard,
    session=session,
    run='01',
    suffix='CT',
    extension='.nii.gz', check=False
).basename

premri_native_bids_fname = BIDSPath(
    subject=subject_wildcard,
    session=session,
    run='01',
    suffix='T1w',
    extension='.nii.gz', check=False
).basename

# output files
# raw T1/CT output
ct_output = os.path.join(BIDS_CT_DIR, ct_native_bids_fname)
t1_output = os.path.join(BIDS_ANAT_DIR, premri_native_bids_fname)

rule prep:
    input:
        # MRI_NIFTI_IMG=expand(t1_output,subject=subjects),
        # CT_bids_fname=expand(ct_output,subject=subjects),
        t1_output=expand(t1_output,subject=subjects),
    params:
        bids_root=bids_root.bids_root,
    log: expand("logs/recon.{subject}.log",subject=subjects)
    container:
        validator_dockerurl
    output:
        report=report('figbids.png',caption='report/figbids.rst',category='BIDS')
    shell:
        "echo 'done';"
        "bids-validator {params.bids_root};"
        "touch figbids.png {output};"

"""
Rule for prepping fs_recon by converting dicoms -> NIFTI images.

For more information, see BIDS specification.
"""

rule convert_dicom_to_nifti_mri:
    params:
        bids_root=bids_root.bids_root,
        subject=subject_wildcard,
        session=session,
        source_folder=sourcedata_path,
    log: "logs/recon.{subject}.log"
    conda:
        "../../envs/heudiconv.yml"
    # container:
    #     heudiconv_dockerurl
    output:
        MRI_bids_fname=t1_output,
    shell:
        "pwd;"\
        "echo {params.session};"
        "heudiconv -d {params.source_folder}/{{subject}}/**/*.dcm -o {params.bids_root} "\
        "-f ./rules/heuristic.py "\
        "-s {params.subject} -ss {params.session} "\
        "-c dcm2niix -b --overwrite --grouping accession_number ;"
