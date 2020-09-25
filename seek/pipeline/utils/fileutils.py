import os
from pathlib import Path

from mne_bids import make_bids_basename
from snakemake.logging import logger

# get the environment variable for freesurfer - for use in getting access to lut's
FREESURFER_HOME = os.getenv("FREESURFER_HOME") or ""
MRTRIX3_HOME = os.getenv("MRTRIX3_HOME") or ""
ATLAS = ["dk", "destrieux"]

# key-word parameters
parc = "aparc.a2009s"
aa = "aparc+aseg"
sval = "pial"
hemispheres = ["lh", "rh"]
resamp_target = "fsaverage5"

BIDS_ROOT = lambda bidsroot: os.getenv("BIDS_ROOT", bidsroot)

DEFAULT_SESSION = "presurgery"
TEMPLATE_SUBJECT = "cvs_avg35_inMNI152"


def _get_subject_center(subjects, centers, subject):
    sub_idx = list(subjects).index(subject)
    return centers[sub_idx]


def _get_session_name(config):
    return config.get("SESSION", DEFAULT_SESSION)


def _get_subject_samples(config):
    import pandas as pd

    samples = pd.read_table(config["subjects"]).set_index("samples", drop=False)


def _get_seek_path():
    """Get relative path to the config file."""
    import seek

    if seek.__file__ is None:
        base_path = None
    else:
        base_path = os.path.dirname(seek.__file__)
    seekdir = os.getenv("SEEKHOME", base_path)

    if seekdir is None:
        raise RuntimeError(
            "Either `seek.__file__` should not be None, "
            "or the environment variable `SEEKHOME` "
            "should be set. Neither is right now."
        )
    logger.debug(f"This is our seek directory: {seekdir}")
    return seekdir


def _get_seek_config():
    config_path = os.path.join(
        _get_seek_path(), "pipeline", "config", "localconfig.yaml"
    )
    return config_path


def ensure_str(func):
    def func_wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        return str(output)

    return func_wrapper


def _get_bids_basename(subject, session, imgtype, ext="nii.gz", **bids_kwargs):
    """Wildcard function to get bids_basename."""
    bids_fname = make_bids_basename(
        subject, session=session, **bids_kwargs, suffix=f"{imgtype}.{ext}"
    )
    return bids_fname


def _get_subject_dir(bids_root, subject):
    return os.path.join(bids_root, f"sub-{subject}")


def _get_anat_bids_dir(bids_root, subject, session):
    return os.path.join(_get_subject_dir(bids_root, subject), f"ses-{session}", "anat")


def _get_ieeg_bids_dir(bids_root, subject, session):
    return os.path.join(_get_subject_dir(bids_root, subject), f"ses-{session}", "ieeg")


def _get_ct_bids_dir(bids_root, subject, session):
    return os.path.join(_get_subject_dir(bids_root, subject), f"ses-{session}", "ct")


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
            self.sourcedir / self.center_id / patient_wildcard / "premri/"
        ).as_posix()

    @ensure_str
    def get_postmri_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.sourcedir / self.center_id / patient_wildcard / "postmri/"
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


SCRIPTS_UTIL_DIR = Path(_get_seek_path()) / "scripts"
