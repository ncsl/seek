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

from seek.pipeline.utils.fileutils import BidsRoot, BIDS_ROOT, _get_session_name, _get_seek_config

configfile: _get_seek_config()

freesurfer_dockerurl = config['freesurfer_docker']
fsl_dockerurl = config['fsl_docker']
blender_dockerurl = config['blender_docker']

BLENDER_PATH = config["blender_path"]
print(BLENDER_PATH)

# get the freesurfer patient directory
bids_root = BidsRoot(BIDS_ROOT(config['bids_root']),
                     center_id=config.get('center_id')
                     )
subject_wildcard = "{subject}"

# initialize directories that we access in this snakemake
FS_DIR = bids_root.freesurfer_dir
FSPATIENT_SUBJECT_FOLDER = bids_root.get_freesurfer_patient_dir(subject_wildcard)
FS_MRI_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "mri"
FS_CT_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "CT"
FS_ELECS_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "elecs"
FS_ACPC_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "acpc"
FS_SURF_FOLDER = Path(FSPATIENT_SUBJECT_FOLDER) / "surf"
FS_LABEL_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER, "label")
FS_ROI_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER, "rois")
FS_OBJ_FOLDER = os.path.join(FSPATIENT_SUBJECT_FOLDER, "obj")

# specific file paths
LH_PIAL_ASC = os.path.join(FS_SURF_FOLDER, "ascii", "lh.pial.asc")
RH_PIAL_ASC = os.path.join(FS_SURF_FOLDER, "ascii", "rh.pial.asc")
LH_ANNOT_FILE = os.path.join(FS_LABEL_FOLDER, "lh.aparc.annot")
RH_ANNOT_FILE = os.path.join(FS_LABEL_FOLDER, "rh.aparc.annot")
LH_PIAL_SRF = os.path.join(FS_SURF_FOLDER, "lh.pial.srf")
RH_PIAL_SRF = os.path.join(FS_SURF_FOLDER, "rh.pial.srf")
LH_ANNOT_DPV = os.path.join(FS_LABEL_FOLDER, "lh.aparc.annot.dpv")
RH_ANNOT_DPV = os.path.join(FS_LABEL_FOLDER, "rh.aparc.annot.dpv")
LH_PIAL_ROI = os.path.join(FS_ROI_FOLDER, "lh.pial_roi")
RH_PIAL_ROI = os.path.join(FS_ROI_FOLDER, "rh.pial_roi")

# blender output file paths
surface_scene_fpath = os.path.join(FSPATIENT_SUBJECT_FOLDER, "blender_objects", "reconstruction.glb")
surface_fbx_fpath = os.path.join(FSPATIENT_SUBJECT_FOLDER, "blender_objects", "reconstruction.fbx")

# coordinate system and electrodes as tsv files
manual_coordsystem_fname = BIDSPath(
    subject=subject_wildcard, session=_get_session_name(config),
    # processing='manual',
    acquisition='seeg', space='fs',
    suffix='coordsystem', extension='.json', root=bids_root.bids_root)
manual_electrodes_fname = BIDSPath(
    subject=subject_wildcard, session=_get_session_name(config),
    # processing='manual',
    acquisition='seeg', space='fs',
    suffix='electrodes', extension='.tsv', root=bids_root.bids_root)

coordsystem_fname = manual_coordsystem_fname #.update(processing='seek')
electrodes_fname = manual_electrodes_fname #.update(processing='seek')

subworkflow contact_labeling_workflow:
    workdir:
           "./rules/"
    snakefile:
             "./rules/label_contacts.smk"
    configfile:
              _get_seek_config()

rule generate_visualization_blender_meshes:
    input:
         surface_scene_file=expand(surface_scene_fpath, subject=subjects),
         # surface_scene_file = os.path.join("./webserver/templates/static/", "reconstruction.glb"),
         # surface_fbx_file = os.path.join("./webserver/templates/static/", "reconstruction.fbx"),
    output:
          report=report('figviz.png', caption='report/figviz.rst', category='Visualization Prep')
    shell:
         "echo 'Done!;'"
         "touch figviz.png {output};"

"""Convert Ascii pial files to surface files."""
rule convert_asc_to_srf:
    input:
         LH_PIAL_ASC=LH_PIAL_ASC,
         RH_PIAL_ASC=RH_PIAL_ASC,
    output:
          LH_PIAL_SRF=LH_PIAL_SRF,
          RH_PIAL_SRF=RH_PIAL_SRF,
    shell:
         "cp {input.LH_PIAL_ASC} {output.LH_PIAL_SRF};"
         "cp {input.RH_PIAL_ASC} {output.RH_PIAL_SRF};"

"""Convert Subcortical surface files to Blender object files."""
rule convert_subcort_to_obj:
    input:
         subcort_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER, "{subject}_subcort_success.txt"),
    params:
          subject=subject_wildcard,
          fsdir=FS_DIR,
          FS_OBJ_FOLDER=str(FS_OBJ_FOLDER),
    container:
        blender_dockerurl
    output:
          obj_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER, "{subject}_subcortobjects_success.txt"),
    shell:
         "export SUBJECTS_DIR={params.fsdir};"
         "mkdir -p {params.FS_OBJ_FOLDER};"
         "./scripts/objMaker.sh {params.subject};"
         "touch {output.obj_success_flag_file};"

"""Convert annotation files to DPV files for Blender."""
rule convert_annot_to_dpv:
    input:
         LH_ANNOT_FILE=LH_ANNOT_FILE,
         RH_ANNOT_FILE=RH_ANNOT_FILE,
    container:
        blender_dockerurl
    output:
          LH_ANNOT_DPV=LH_ANNOT_DPV,
          RH_ANNOT_DPV=RH_ANNOT_DPV,
    shell:
         "./scripts/annot2dpv {input.LH_ANNOT_FILE} {output.LH_ANNOT_DPV};"
         "./scripts/annot2dpv {input.RH_ANNOT_FILE} {output.RH_ANNOT_DPV};"

"""Split surfaces into files per right/left hemisphere."""
rule split_surfaces:
    input:
         LH_ANNOT_DPV=LH_ANNOT_DPV,
         RH_ANNOT_DPV=RH_ANNOT_DPV,
         LH_PIAL_SRF=LH_PIAL_SRF,
         RH_PIAL_SRF=RH_PIAL_SRF,
    params:
          LH_PIAL_ROI=LH_PIAL_ROI,
          RH_PIAL_ROI=RH_PIAL_ROI,
    container:
        blender_dockerurl
    output:
          roi_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER, "surfaces_roi_flag_success.txt")
    shell:
         # "cd ./scripts/;"
         "./scripts/splitsrf {input.LH_PIAL_SRF} {input.LH_ANNOT_DPV} {params.LH_PIAL_ROI};"
         "./scripts/splitsrf {input.RH_PIAL_SRF} {input.RH_ANNOT_DPV} {params.RH_PIAL_ROI};"
         "touch {output.roi_flag_file};"

rule create_surface_objects:
    input:
         roi_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER, "surfaces_roi_flag_success.txt"),
         electrode_fpath=contact_labeling_workflow(electrodes_fname),
         obj_success_flag_file=os.path.join(FSPATIENT_SUBJECT_FOLDER, "{subject}_subcortobjects_success.txt"),
    params:
          LH_PIAL_ROI=LH_PIAL_ROI,
          RH_PIAL_ROI=RH_PIAL_ROI,
          fsdir=FS_DIR,
          subject=subject_wildcard,
          materialcolors_file=os.path.join(os.getcwd(), "./scripts/materialColors.json"),
          BLENDER_PATH=BLENDER_PATH,
    container:
        blender_dockerurl
    output:
          surface_scene_file=surface_scene_fpath,
          surface_fbx_file=surface_fbx_fpath,
    shell:
         "echo 'Creating surface objects for rendering!';"
         "export SUBJECTS_DIR={params.fsdir};"
         "./scripts/surfaceToObject.sh {params.subject};"
         "{params.BLENDER_PATH} --background --python ./scripts/sceneCreator.py -- " \
                             "{params.fsdir} " \
                             "{params.subject} " \
                             "{input.electrode_fpath} " \
                             "True False " \
                             "{output.surface_fbx_file} " \
                             "{output.surface_scene_file} " \
                             "{params.materialcolors_file};"
         # "cp {output.surface_scene_file} {output.copied_surface_scene_file};"
         # "cp {output.surface_fbx_file} {output.copied_surface_fbx_file};"
