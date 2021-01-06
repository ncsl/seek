import os
from pathlib import Path

from mne_bids import BIDSPath
from snakemake.logging import logger

# get the environment variable for freesurfer - for use in getting access to lut's
FREESURFER_HOME = os.getenv("FREESURFER_HOME") or ""
MRTRIX3_HOME = os.getenv("MRTRIX3_HOME") or ""
SCRIPTS_UTIL_DIR = "../format"
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
    logger.debug(f"Found {centers[sub_idx]} for {subject}")
    return centers[sub_idx]


def _get_session_name(config):
    return config.get("SESSION", DEFAULT_SESSION)


def _get_subject_samples(config):
    import pandas as pd

    samples = pd.read_table(config["subjects"]).set_index("samples", drop=False)


def _get_seek_path():
    """Get relative path to the config file."""
    import seek

    logger.debug(f"Attempting to infer seek filepath: {seek.__file__}.")
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
    config_path = Path(_get_seek_path()).rglob("config/localconfig.yaml")
    config_path = list(config_path)[0]
    logger.info(f"Found configuration filepath: {config_path}")
    # config_path = os.path.join(
    #     _get_seek_path(), 'seek', "pipeline", "config", "localconfig.yaml"
    # )
    return config_path


def ensure_str(func):
    def func_wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        return str(output)

    return func_wrapper


def _get_bids_basename(subject, session, imgtype, ext="nii.gz", **bids_kwargs):
    """Wildcard function to get bids_basename."""
    bids_fname = BIDSPath(
        subject,
        session=session,
        **bids_kwargs,
        suffix=imgtype,
        extension=ext,
        check=False,
    ).basename
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

    def __init__(
        self,
        subject_id,
        bids_root,
        site_id=None,
        derivatives_dir=None,
        subject_wildcard="{subject}",
        **kwargs,
    ):
        self.subject_id = subject_id

        bids_root = Path(bids_root)
        if not bids_root.exists():
            raise RuntimeError(f"Bidsroot {bids_root} does not exist.")
        self.bids_root = bids_root
        self.site_id = site_id
        self.subject_wildcard = subject_wildcard

        # derivatives directory is either custom, or derived from bids_root.
        if derivatives_dir is None:
            self.derivatives_dir = Path(self.bids_root / "derivatives")
        else:
            self.derivatives_dir = Path(derivatives_dir)

    def __repr__(self):
        return self.bids_root

    @property
    def freesurfer_dir(self):
        return Path(self.derivatives_dir / "freesurfer")

    @property
    def elecs_dir(self):
        return Path(self.freesurfer_dir / self.subject_id / "elecs")

    @property
    def mesh_dir(self):
        return Path(self.freesurfer_dir / self.subject_id / "Meshes")

    @property
    def gyri_dir(self):
        return Path(self.freesurfer_dir / self.subject_id / "label" / "gyri")

    @property
    def mri_dir(self):
        return Path(self.freesurfer_dir / self.subject_id / "mri")

    @property
    def source_chain(self):
        return self._get_source_chain()

    def _get_source_chain(self):
        source_chain = Path("sourcedata")
        if self.site_id is not None:
            source_chain = source_chain / self.site_id
        source_chain = source_chain / self.subject_wildcard
        return source_chain

    def get_surface_atlas_suffix(self, atlas_type):
        if atlas_type == "desikan-killiany":
            surf_atlas_suffix = ""
        elif atlas_type == "destrieux":
            surf_atlas_suffix = "--a2009s"
        else:
            surf_atlas_suffix = ""
        return surf_atlas_suffix

    def get_depth_atlas_suffix(self, atlas_type):
        if atlas_type == "desikan-killiany":
            depth_atlas_suffix = ""
        elif atlas_type == "destrieux":
            depth_atlas_suffix = ".a2009s"
        else:
            depth_atlas_suffix = ".a2009s"
        return depth_atlas_suffix

    def get_freesurfer_patient_dir(self):
        return Path(self.derivatives_dir / "freesurfer" / self.subject_wildcard)

    def get_premri_dir(self):
        return Path(self.bids_root / self.source_chain / "premri").as_posix()

    def get_postmri_dir(self):
        return Path(self.bids_root / self.source_chain / "postmri").as_posix()

    def get_rawct_dir(self):
        return Path(self.bids_root / self.source_chain / "postct").as_posix()

    def get_rawacpc_dir(self):
        return Path(self.bids_root / self.source_chain / "acpc").as_posix()
