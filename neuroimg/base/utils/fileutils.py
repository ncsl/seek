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


class PatientBidsRoot:
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

    def __init__(self, subject_id, bids_root, derivatives_dir=None):
        self.subject_id = subject_id

        bids_root = Path(bids_root)
        if not bids_root.exists():
            raise RuntimeError(f"Bidsroot {bids_root} does not exist.")
        self.bids_root = bids_root

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

    def get_freesurfer_patient_dir(self, patient_wildcard="{subject}"):
        return Path(self.derivatives_dir / "freesurfer" / patient_wildcard)

    def get_premri_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.bids_root / "sourcedata" / "neuroimaging" / patient_wildcard / "premri"
        ).as_posix()

    def get_postmri_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.bids_root
            / "sourcedata"
            / "neuroimaging"
            / patient_wildcard
            / "postmri"
        ).as_posix()

    def get_rawct_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.bids_root / "sourcedata" / "neuroimaging" / patient_wildcard / "postct"
        ).as_posix()

    def get_rawacpc_dir(self, patient_wildcard="{subject}"):
        return Path(
            self.bids_root / "sourcedata" / "neuroimaging" / patient_wildcard / "acpc"
        ).as_posix()
