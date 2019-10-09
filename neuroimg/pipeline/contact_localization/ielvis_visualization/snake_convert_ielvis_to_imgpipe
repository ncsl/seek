import os


"""
snakemake --dag | dot -Tpdf > dag_neuroimaging_pipeline_contactlocalization.pdf
"""

configfile: "../config/localconfig.yaml"

# get the environment variable for freesurfer - for use in getting access to lut's
HOME_DIR="../../"
SCRIPTS_UTIL_DIR = "../../neuroimg/format/"
ATLAS=['dk', 'destrieux']

''' USER DEFINED DIRECTORIES TO STORE FINAL DATA IN ORGANIZED SUBFOLDERS '''
FS_PATIENT_OUTPUT_DIR = os.path.join(config['fs_outputdatadir'],
                                     "{patient_id}")

COREGISTRATION_OUTPUT_DIR = os.path.join(config['fs_outputdatadir'],
                                         "{patient_id}",
                                         "coregistration")

ELECSDIR = os.path.join(FS_PATIENT_OUTPUT_DIR, "elecs")

IELVIS_OUTPUT_DIR = os.path.join(FS_PATIENT_OUTPUT_DIR, "elec_recon")


# First rule
rule all:
    input:
        elecfile = expand(os.path.join(IELVIS_OUTPUT_DIR,
                                                      "brainmask.mgz"),
                                                patient_id = config['patients']),
    shell:
        "echo 'done'"



rule add_ras_data:
    input:
        MRI_NIFTI_IMG = os.path.join(FS_PATIENT_OUTPUT_DIR,
                             "mri",
                             "T1.nii"),
        eleccoords_file = os.path.join(ELECSDIR,
                                       "{patient_id}_elecs_all.mat"),

    output:
        outputfile = os.path.join(ELECSDIR,
                               "{patient_id}_elecs_all_RAS.mat"),
    shell:
        "python -m convert_ielvis_to_imgpipe {input.MRI_NIFTI_IMG} " \
                                        "{input.eleccoords_file} " \
                                        "{output.outputfile};"

"""
"""
rule move_files_to_ielvis:
    input:
        elecfile = os.path.join(ELECSDIR,
                          "{patient_id}_elecs_all_RAS.mat"),
        # DEPENDENCY ON RECONSTRUCTION WORKFLOW
        CT_NIFTI_IMG = os.path.join(FS_PATIENT_OUTPUT_DIR, "CT",
                                                       "CT.nii.gz"),
        brainmask_file = os.path.join(FS_PATIENT_OUTPUT_DIR ,
                                                                     "mri",
                                                                     "brainmask.mgz"),
        # mapped image from CT -> MRI
        CT_IN_PRE_NIFTI_IMG = os.path.join(FS_PATIENT_OUTPUT_DIR, "CT",
                                                      "rCT.nii.gz"),
        MRI_NIFTI_IMG = os.path.join(FS_PATIENT_OUTPUT_DIR,
                                     "mri",
                                                    "T1.nii"),
    params:
        IELVISDIR = IELVIS_OUTPUT_DIR,
    output:
        CT_IN_PRE_NIFTI_IMG = os.path.join(IELVIS_OUTPUT_DIR,
                                                           "CT_IN_T1.nii.gz"),
        MRI_NIFTI_IMG = os.path.join(IELVIS_OUTPUT_DIR,
                             "T1.nii"),
        brainmask_file = os.path.join(IELVIS_OUTPUT_DIR,
                                                      "brainmask.mgz"),
        elecfile = os.path.join(IELVIS_OUTPUT_DIR,
                        "{patient_id}_elecs_all_RAS.mat"),
    shell:
        # "mkdir {params.IELVISDIR};"
        "cp {input.CT_IN_PRE_NIFTI_IMG} {output.CT_IN_PRE_NIFTI_IMG};"
        "cp {input.MRI_NIFTI_IMG} {output.MRI_NIFTI_IMG};"
        "cp {input.brainmask_file} {output.brainmask_file};"
        "cp {input.elecfile} {output.elecfile};"

