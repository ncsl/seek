import os

# get the environment variable for freesurfer - for use in getting access to lut's
FREESURFER_HOME = os.getenv("FREESURFER_HOME") or ''
MRTRIX3_HOME = os.getenv("MRTRIX3_HOME") or ''
SCRIPTS_UTIL_DIR = "../../format"
ATLAS=['dk', 'destrieux']

# key-word parameters
parc = "aparc.a2009s"
aa = "aparc+aseg"
sval = "pial"
hemispheres = ['lh', 'rh']
resamp_target = "fsaverage5"

def get_freesurfer_dir(config):
    """
    Function to return the FreeSurfer directory based on a config file. Wraps a wildcard for patient_id.

    Parameters
    ----------
    config :
    patient_wildcard :

    Returns
    -------

    """
    ''' USER DEFINED DIRECTORIES TO STORE FINAL DATA IN ORGANIZED SUBFOLDERS '''
    FS_PATIENT_OUTPUT_DIR = os.path.join(config['basedatadir'],
                                         "derivatives",
                                         "freesurfer")
    return FS_PATIENT_OUTPUT_DIR

def get_freesurfer_patient_dir(config, patient_wildcard="{patient_id}"):
    """
    Function to return the FreeSurfer directory based on a config file. Wraps a wildcard for patient_id.

    Parameters
    ----------
    config :
    patient_wildcard :

    Returns
    -------

    """
    ''' USER DEFINED DIRECTORIES TO STORE FINAL DATA IN ORGANIZED SUBFOLDERS '''
    FS_PATIENT_OUTPUT_DIR = os.path.join(config['basedatadir'],
                                         "derivatives",
                                         "freesurfer",
                                         patient_wildcard)
    return FS_PATIENT_OUTPUT_DIR

def get_rawmri_dir(config, patient_wildcard="{patient_id}"):
    return os.path.join(config['rawdatadir'],
                 patient_wildcard,
                 "premri")

def get_rawct_dir(config, patient_wildcard="{patient_id}"):
    return os.path.join(config['rawdatadir'],
                 patient_wildcard,
                 "postct")

def get_rawacpc_dir(config, patient_wildcard="{patient_id}"):
    return os.path.join(config['rawdatadir'],
                 patient_wildcard,
                 "acpc")