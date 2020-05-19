import os
from pathlib import Path

# get the environment variable for freesurfer - for use in getting access to lut's
FREESURFER_HOME = os.getenv("FREESURFER_HOME") or ""
MRTRIX3_HOME = os.getenv("MRTRIX3_HOME") or ""
SCRIPTS_UTIL_DIR = "../../format"
ATLAS = ["dk", "destrieux"]

# key-word parameters
parc = "aparc.a2009s"
aa = "aparc+aseg"
sval = "pial"
hemispheres = ["lh", "rh"]
resamp_target = "fsaverage5"

BIDS_ROOT = lambda bidsroot: os.getenv("BIDS_ROOT", bidsroot)

DEFAULT_SESSION = "presurgery"

def _get_session_name(config):
    return config.get('session', DEFAULT_SESSION)

def _get_seek_config():
    """Get relative path to the config file."""
    import seek

    base_path = os.path.dirname(seek.__file__)
    seekdir = os.getenv("SEEKHOME", base_path)
    config_path = os.path.join(seekdir, "pipeline", "config", "localconfig.yaml")
    return config_path


def ensure_str(func):
    def func_wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        return str(output)

    return func_wrapper


class BidsRoot:
    """
    A utility class wrapper for a Bids root.

    Bids root inherits the BIDS standard for neuroimaging data, and
    stores data within the bids-root in a very strict standard. See https://bids.neuroimaging.io/
    for more information on BIDS.

    BidsRoot stores FreeSurfer (T1 MRI) derived data within "derivatives/" sub-folder.
    It also contains specifically laid out sourcedata in the form of "dicoms" for
    pre-implantation T1 MRI, post-implantation CT and post-surgery T1 MRI.

    premri is used to run FreeSurfer on the patient brain.
    postct is used to localize implanted electrode contacts.
    postmri is used to estimate the region of surgical intervention.

    """

    def __init__(self, bids_root, derivatives_dir=None, center_id=None):
        bids_root = Path(bids_root)
        if not bids_root.exists():
            raise RuntimeError(f"Bidsroot {bids_root} does not exist.")
        self.bids_root = bids_root
        self.sourcedir = Path(self.bids_root / "sourcedata")

        if center_id is None:
            center_id = ""
        self.center_id = center_id

        # derivatives directory is either custom, or derived from bids_root.
        if derivatives_dir is None:
            self.derivatives_dir = Path(self.bids_root / "derivatives")
        else:
            self.derivatives_dir = Path(derivatives_dir)

    def __repr__(self):
        return self.bids_root

    # def _check_center_id(self, patient_wildcard):
    #     source_children = [x for x in self.sourcedir.glob("*") if patient_wildcard

    @property
    def freesurfer_dir(self):
        return Path(self.derivatives_dir / "freesurfer")

    @ensure_str
    def get_freesurfer_patient_dir(self, patient_wildcard="{subject}"):
        return Path(self.derivatives_dir / "freesurfer" / patient_wildcard).as_posix()

    @ensure_str
    def get_premri_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.sourcedir / self.center_id / patient_wildcard / "premri"
        ).as_posix()

    @ensure_str
    def get_postmri_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.sourcedir / self.center_id / patient_wildcard / "postmri"
        ).as_posix()

    @ensure_str
    def get_rawct_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.sourcedir / self.center_id / patient_wildcard / "postct"
        ).as_posix()

    @ensure_str
    def get_rawacpc_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.sourcedir / self.center_id / patient_wildcard / "acpc"
        ).as_posix()
