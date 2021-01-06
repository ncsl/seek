"""
===============================================
05. Visualization workflow on web browser
================================================

This pipeline depends on the following functions:

    * blender
    * Flask

from Blender, Flask
"""

import os
import sys
from pathlib import Path

from mne_bids import BIDSPath

sys.path.append("../../../")

from seek.utils.fileutils import (BidsRoot, BIDS_ROOT, _get_seek_path,
                                  _get_session_name, _get_seek_config)

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']
fsl_dockerurl = config['fsl_docker']
blender_dockerurl = config['blender_docker']

# get the freesurfer patient directory
subject_wildcard = "{subject}"
bids_root = BidsRoot(subject_wildcard,BIDS_ROOT(config['bids_root']),
    site_id=config['site_id'],subject_wildcard=subject_wildcard)

SEEKHOME = _get_seek_path()
scripts_dir = os.path.join(SEEKHOME,'workflow','prep_vizengine_workflow','scripts')

# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_FOLDER = bids_root.get_freesurfer_patient_dir()
FS_MRI_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "mri"
FS_CT_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "CT"
FS_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "elecs"
FS_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "acpc"
FS_SURF_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "surf"
FS_LABEL_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER,"label")
FS_ROI_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER,"rois")
FS_OBJ_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER,"obj")

# specific file paths
LH_PIAL_ASC = os.path.join(FS_SURF_FOLDER,"ascii","lh.pial.asc")
RH_PIAL_ASC = os.path.join(FS_SURF_FOLDER,"ascii","rh.pial.asc")
LH_ANNOT_FILE = os.path.join(FS_LABEL_FOLDER,"lh.aparc.annot")
RH_ANNOT_FILE = os.path.join(FS_LABEL_FOLDER,"rh.aparc.annot")
LH_PIAL_SRF = os.path.join(FS_SURF_FOLDER,"lh.pial.srf")
RH_PIAL_SRF = os.path.join(FS_SURF_FOLDER,"rh.pial.srf")
LH_ANNOT_DPV = os.path.join(FS_LABEL_FOLDER,"lh.aparc.annot.dpv")
RH_ANNOT_DPV = os.path.join(FS_LABEL_FOLDER,"rh.aparc.annot.dpv")
LH_PIAL_ROI = os.path.join(FS_ROI_FOLDER,"lh.pial_roi")
RH_PIAL_ROI = os.path.join(FS_ROI_FOLDER,"rh.pial_roi")

# blender output file paths
surface_scene_fpath = os.path.join(FSPATIENT_SUBJECT_FOLDER,"blender_objects", "brain.glb")
surface_fbx_fpath = os.path.join(FSPATIENT_SUBJECT_FOLDER,"blender_objects","brain.fbx")
electrodes_scene_fpath = os.path.join(FSPATIENT_SUBJECT_FOLDER,"blender_objects","electrodes.glb")

# coordinate system and electrodes as tsv files
manual_coordsystem_fname = BIDSPath(
    subject=subject_wildcard,session=_get_session_name(config),
    # processing='manual',
    acquisition='seeg',space='fs',datatype='ieeg',
    suffix='coordsystem',extension='.json',root=bids_root.bids_root)
manual_electrodes_fname = BIDSPath(
    subject=subject_wildcard,session=_get_session_name(config),
    # processing='manual',
    datatype='ieeg',
    acquisition='seeg',space='fs',
    suffix='electrodes',extension='.tsv',root=bids_root.bids_root)

coordsystem_fname = manual_coordsystem_fname  #.update(processing='seek')
electrodes_fname = manual_electrodes_fname  #.update(processing='seek')

subworkflow contact_labeling_workflow:
    workdir:
        "../contact_labeling_workflow/"
    snakefile:
        "../contact_labeling_workflow/Snakefile"
    configfile:
        _get_seek_config()

rule generate_visualization_blender_meshes:
    input:
        surface_scene_file=expand(surface_scene_fpath,subject=subjects),
    # electrodes_scene_file=expand(electrodes_scene_fpath,subject=subjects),
    # surface_scene_file = os.path.join("./webserver/templates/static/", "reconstruction.glb"),
    # surface_fbx_file = os.path.join("./webserver/templates/static/", "reconstruction.fbx"),
    log:
        expand("logs/prep_vizengine.{subject}.log",subject=subjects)
    output:
        report=report('figviz.png',caption='report/figviz.rst',category='Visualization Prep')
    shell:
        "echo 'Done!';"
        "touch figviz.png {output};"

"""Convert Ascii pial files to surface files."""

rule convert_asc_to_srf:
    input:
        LH_PIAL_ASC=LH_PIAL_ASC,
        RH_PIAL_ASC=RH_PIAL_ASC,
    log:
        "logs/prep_vizengine.{subject}.log"
    output:
        LH_PIAL_SRF=LH_PIAL_SRF,
        RH_PIAL_SRF=RH_PIAL_SRF,
    shell:
        "cp {input.LH_PIAL_ASC} {output.LH_PIAL_SRF};"
        "cp {input.RH_PIAL_ASC} {output.RH_PIAL_SRF};"

"""Convert annotation files to DPV files for Blender."""

rule convert_annot_to_dpv:
    input:
        LH_ANNOT_FILE=LH_ANNOT_FILE,
        RH_ANNOT_FILE=RH_ANNOT_FILE,
    params:
        scripts_dir=scripts_dir
    log:
        "logs/prep_vizengine.{subject}.log"
    container:
        blender_dockerurl
    output:
        LH_ANNOT_DPV=LH_ANNOT_DPV,
        RH_ANNOT_DPV=RH_ANNOT_DPV,
    shell:
        "cd {params.scripts_dir}/octave/;"  # needed to have `read_annotation.m` in path
        "{params.scripts_dir}/octave/annot2dpv {input.LH_ANNOT_FILE} {output.LH_ANNOT_DPV};"
        "{params.scripts_dir}/octave/annot2dpv {input.RH_ANNOT_FILE} {output.RH_ANNOT_DPV};"

"""Split surfaces into files per right/left hemisphere."""

rule split_surfaces:
    input:
        LH_ANNOT_DPV=LH_ANNOT_DPV,
        RH_ANNOT_DPV=RH_ANNOT_DPV,
        LH_PIAL_SRF=LH_PIAL_SRF,
        RH_PIAL_SRF=RH_PIAL_SRF,
    log:
        "logs/prep_vizengine.{subject}.log"
    params:
        LH_PIAL_ROI=LH_PIAL_ROI,
        RH_PIAL_ROI=RH_PIAL_ROI,
        scripts_dir=scripts_dir,
        FS_ROI_FOLDER=FS_ROI_FOLDER,
    container:
        blender_dockerurl
    output:
        roi_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"surfaces_roi_flag_success.txt")
    shell:
        "mkdir -p {params.FS_ROI_FOLDER};"
        "cd {params.scripts_dir}/octave/;"  # needed to have `srfread/srfwrite` in path
        "echo '{params.scripts_dir}/octave/splitsrf {input.LH_PIAL_SRF} {input.LH_ANNOT_DPV} {params.LH_PIAL_ROI}';"
        "{params.scripts_dir}/octave/splitsrf {input.LH_PIAL_SRF} {input.LH_ANNOT_DPV} {params.LH_PIAL_ROI};"
        "{params.scripts_dir}/octave/splitsrf {input.RH_PIAL_SRF} {input.RH_ANNOT_DPV} {params.RH_PIAL_ROI};"
        "touch {output.roi_flag_file};"

"""Convert Subcortical surface files to Blender object files."""

rule convert_subcort_to_blenderobj:
    input:
        subcort_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"{subject}_subcort_success.txt"),
    params:
        subject=subject_wildcard,
        fsdir=FS_DIR,
        FS_OBJ_FOLDER=str(FS_OBJ_FOLDER),
        scripts_dir=scripts_dir
    log:
        "logs/prep_vizengine.{subject}.log"
    container:
        blender_dockerurl
    output:
        obj_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"{subject}_subcortobjects_success.txt"),
    shell:
        "export SUBJECTS_DIR={params.fsdir};"
        "mkdir -p {params.FS_OBJ_FOLDER};"
        "{params.scripts_dir}/bash/objMaker.sh {params.subject};"
        "touch {output.obj_success_flag_file};"

"""Rule to create blender ``.obj`` files from ``.srf`` files."""

rule convert_cortical_to_blenderobj:
    input:
        roi_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"surfaces_roi_flag_success.txt"),
        obj_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"{subject}_subcortobjects_success.txt"),
    params:
        LH_PIAL_ROI=LH_PIAL_ROI,
        RH_PIAL_ROI=RH_PIAL_ROI,
        subject=subject_wildcard,
        scripts_dir=scripts_dir,
        fsdir=FS_DIR,
    log:
        "logs/prep_vizengine.{subject}.log"
    container:
        blender_dockerurl
    output:
        surface_obj_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"surfaces_obj_flag_success.txt"),
    shell:
        "echo 'Creating surface objects for rendering!';"
        "export SUBJECTS_DIR={params.fsdir};"
        "{params.scripts_dir}/bash/surfaceToObject.sh {params.subject};"
        "touch {output.surface_obj_flag_file};"


"""Rule to create brain surface ``.glb`` files."""

rule create_brain_glb_files:
    input:
        surface_obj_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"surfaces_obj_flag_success.txt"),
        obj_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER,"{subject}_subcortobjects_success.txt"),
    params:
        LH_PIAL_ROI=LH_PIAL_ROI,
        RH_PIAL_ROI=RH_PIAL_ROI,
        fsdir=FS_DIR,
        subject=subject_wildcard,
        materialcolors_file=os.path.join(scripts_dir, "octave/materialColors.json"),
        scripts_dir=scripts_dir,
    log:
        "logs/prep_vizengine.{subject}.log"
    container:
        blender_dockerurl
    output:
        surface_scene_file=surface_scene_fpath,
        # surface_fbx_file=surface_fbx_fpath,
    shell:
        "echo 'Creating brain glb objects for rendering!';"
        "export SUBJECTS_DIR={params.fsdir};"
        "export SUBJECT={params.subject};"
        "/usr/local/blender/blender --background {params.scripts_dir}/startup.blend " \
        "--python {params.scripts_dir}/brain_generator.py " \
        "-- {params.materialcolors_file};"


"""Rule to create electrode in brain coordinate system ``.glb`` files."""

rule create_electrode_glb_files:
    input:
        electrode_fpath=str(electrodes_fname)
    params:
        fsdir=FS_DIR,
        subject=subject_wildcard,
        materialcolors_file=os.path.join(scripts_dir, "octave/materialColors.json"),
        scripts_dir=scripts_dir
    log:
        "logs/prep_vizengine.{subject}.log"
    container:
        blender_dockerurl
    output:
        electrodes_scene_file=electrodes_scene_fpath,
    shell:
        "echo 'Creating brain glb objects for rendering!';"
        "export SUBJECTS_DIR={params.fsdir};"
        "export SUBJECT={params.subject};"
        "blender --background {params.scripts_dir}/startup.blend " \
        "--python {params.scripts_dir}/electrode_generator.py " \
        "-- {input.electrode_fpath};"
