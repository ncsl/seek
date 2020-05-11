import collections
import logging
import re
from typing import Dict, List, Tuple, Union

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from natsort import natsorted
from nibabel.affines import apply_affine

logger = logging.getLogger(__name__)


class Contact:
    """
    Contact object, which is the point of an iEEG recording.

    Parameters
    ----------
    ch_name : str
    coord : tuple
    coord_type : str
    coordsystem : str

    """

    def __init__(
        self,
        ch_name: str,
        coord: Union[Tuple, List, np.ndarray],
        coord_type: str,
        coordsystem: str = None,
    ):
        self.name = ch_name

        electrode, ch_num = re.match("^([A-Za-z]+[']?)([0-9]+)$", self.name).groups()
        self.electrode = electrode
        self.number = int(ch_num)
        self.coord = coord
        self.coordsystem = coordsystem
        self.coord_type = coord_type

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def set_coord_system(self, coordsystem: str):
        """Set coordinate system of electrodes (e.g. RAS)."""
        self.coordsystem = coordsystem

    def transform_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Transform coordinates between mm and voxel space."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )

        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrode {self.name} is "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )
        new_coord = self.get_transformed_coords(img, coord_type)
        self.coord = new_coord
        self.coord_type = coord_type

    def get_transformed_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Get transformed coordinates of a certain type."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )

        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrode {self.name} is "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )

        affine = img.affine

        # assign the affine transformation to use depending
        # if going from vox -> mm, or mm -> vox
        if coord_type == "vox":
            inv_affine = npl.inv(affine)
        else:
            inv_affine = affine

        vox = apply_affine(inv_affine, self.coord)
        return vox


class Electrode:
    """
    Electrode object consisting of many contacts along the electrode.

    This is primarily for SEEG, but could be used for ECoG as well.
    Each electrode might have anywhere from 6 - 20 contacts on it.

    Parameters
    ----------
    ch_names :
    ch_coords :
    coord_type :
    """

    def __init__(self, ch_names: List, ch_coords: List, coord_type: str):
        if len(ch_coords) != len(ch_names):
            raise RuntimeError(
                "Need to pass in the same amount of coordinates as channels. "
                f"Passed in {len(ch_names)} chs and {len(ch_coords)} coords."
            )

        self.contacts = natsorted(
            [
                Contact(ch_name, ch_coord, coord_type)
                for ch_name, ch_coord in zip(ch_names, ch_coords)
            ],
            key=str,
        )

        _elecs = set([x.electrode for x in self.contacts])
        if len(_elecs) > 1:
            raise RuntimeError(
                f"Number of electrodes in passed data is {len(_elecs)}. "
                f"Each Electrode should only have."
            )

        self.name = self.contacts[0].electrode

    def __repr__(self):
        return self.name

    def __str__(self):
        return ",".join([x.name for x in self.contacts])

    def __len__(self):
        return len(self.contacts)

    @property
    def coord_type(self):
        """Coordinate type (mm, voxel) that electrode coordinates are in."""
        coord_type = np.unique([x.coord_type for x in self.contacts])
        if len(coord_type) > 1:
            raise RuntimeError(
                "Inside Electrode, " "Coordinate type should not be more then 1..."
            )
        return coord_type[0]

    def remove_contact(self, ch_name):
        """Remove contact from electrode."""
        contacts = []
        for contact in self.contacts:
            if contact.name != ch_name:
                contacts.append(contact)
            else:
                logger.info(f"Removed contact {ch_name}!")
        self.contacts = natsorted(contacts, key=str)

    def get_entry_ch(self):
        """Get the initial channel (most medial)."""
        return self.contacts[0]

    def get_exit_ch(self):
        """Get the channel at the end of the electrode."""
        return self.contacts[-1]

    def transform_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Transform coordinates."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )
        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrode {self.name} is "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )
        for contact in self.contacts:
            contact.transform_coords(img, coord_type)

    def get_transformed_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Get transformed coordinates."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )
        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrode {self.name} is "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )

        voxels = []
        for contact in self.contacts:
            voxels.append(contact.get_transformed_coords(img, coord_type))
        return voxels


class ElectrodeIterator:
    """Iterator class for Electrode."""

    def __init__(self, electrodes):
        # Team object reference
        self._electrodes = electrodes
        # member variable to keep track of current index
        self._index = 0

    def __next__(self):
        """Get the next value from team object's lists."""
        if self._index < (len(self._electrodes)):
            result = self._electrodes[self._index]
            self._index += 1
            return result
        # End of Iteration
        raise StopIteration


class Electrodes:
    """
    Class object for a set of Electrodes.

    Electrodes consist of many Electrode objects.

    Parameters
    ----------
    ch_names :
    ch_coords :
    coord_type :
    """

    def __init__(self, ch_names: List, ch_coords: List, coord_type: str):
        if len(ch_coords) != len(ch_names):
            raise RuntimeError(
                "Need to pass in the same amount of coordinates as channels. "
                f"Passed in {len(ch_names)} chs and {len(ch_coords)} coords."
            )

        self._init_elecs(ch_names, ch_coords, coord_type)

    def __repr__(self):
        return self.electrodes

    def __str__(self):
        return ",".join([x.name for x in self.electrodes])

    def __iter__(self):
        return ElectrodeIterator(self.electrodes)

    def __len__(self):
        return len(self.electrodes)

    def _init_elecs(self, ch_names, ch_coords, coord_type):
        """Initialize class."""
        electrodes = []

        # create a dictionary of electrodes -> contacts -> coordinates
        labeled_electrodes = {}
        for name, coord in zip(ch_names, ch_coords):
            elecname = re.findall(r"[A-Za-z']+", name)[0]
            labeled_electrodes.setdefault(elecname, {})[name] = coord

        #
        for elecname, elec in labeled_electrodes.items():
            ch_names = list(elec.keys())
            ch_coords = list(elec.values())
            electrode = Electrode(ch_names, ch_coords, coord_type)
            electrodes.append(electrode)

        self.electrodes = electrodes

    @property
    def coord_type(self):
        """Get coordinate type."""
        coord_type = np.unique([x.coord_type for x in self.electrodes])
        if len(coord_type) > 1:
            raise RuntimeError(
                "Inside Electrodes, " "Coordinate type should not be more then 1..."
            )
        return coord_type[0]

    def get_electrode(self, elecname):
        """Get an electrode by its name."""
        for elec in self.electrodes:
            if elec.name == elecname:
                return elec

        raise RuntimeError(f"No electrode {elecname} present.")

    def get_voxels(self, img: nb.Nifti2Image) -> Dict:
        """Get voxels of the electrode coordinates."""
        voxel_dict = collections.defaultdict(list)
        for electrode in self.electrodes:
            voxel_dict[electrode.name].extend(electrode.get_transformed_coords(img))
        return voxel_dict

    def get_transformed_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Get transformed coordinates."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )
        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrodes are "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )

        voxel_dict = collections.defaultdict(list)
        for electrode in self.electrodes:
            voxel_dict[electrode.name].extend(
                electrode.get_transformed_coords(img, coord_type)
            )
        return voxel_dict

    def transform_coords(self, img: nb.Nifti2Image, coord_type: str):
        """Transform coordinates into mm, or voxel space."""
        if coord_type not in ["mm", "vox"]:
            raise ValueError(
                "Accepted coord_type are 'mm', or 'vox'. "
                f"You passed in {coord_type}."
            )
        if self.coord_type == coord_type:
            raise RuntimeError(
                f"Currently Electrodes are "
                f"already in {self.coord_type} space "
                f"Pass in `coord_type` = 'mm', or 'vox'."
            )

        for electrode in self.electrodes:
            electrode.transform_coords(img, coord_type)
