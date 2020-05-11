import enum
import os

import nibabel
import numpy as np

"""
Stores the base classes for neuroimaging interaction.

- Hemisphere: an enum object to denote which side of the brain 
- ColorLut: wrapper to access the color LUT.
- RegionIndexMapping: wrapper to map region indices from one LUT -> another LUT

"""

SUBCORTICAL_REG_INDS = [
    8,
    10,
    11,
    12,
    13,
    16,
    17,
    18,
    26,
    47,
    49,
    50,
    51,
    52,
    53,
    54,
    58,
]
FS_LUT_LH_SHIFT = 1000
FS_LUT_RH_SHIFT = 2000


class Hemisphere(enum.Enum):
    """
    An enumeration of the brain hemisphere.

    Only includes right and left side.

    """

    rh = "rh"
    lh = "lh"


class RegionIndexMapping(object):
    """
    Class wrapper for a region index mapping.

    This maps each index of the source file to an index in the corresponding look up table.

    Parameters
    ----------
    color_lut_src_file: The source lookup table
    color_lut_trg_file: The target's lookup table.

    """

    def __init__(
        self, color_lut_src_file: os.PathLike, color_lut_trg_file: os.PathLike
    ):
        self.src_table = ColorLut(color_lut_src_file)
        self.trg_table = ColorLut(color_lut_trg_file)

        names_to_trg = dict(zip(self.trg_table.names, self.trg_table.inds))

        self.src_to_trg = dict()
        for src_ind, src_name in zip(self.src_table.inds, self.src_table.names):
            trg_ind = names_to_trg.get(src_name, None)
            if trg_ind is not None:
                self.src_to_trg[src_ind] = trg_ind

        self.unknown_ind = names_to_trg.get(
            "Unknown", 0
        )  # zero as the default unknown area

    def source_to_target(self, index):
        """
        Convert source index on the Source LUT to the Target LUT index.

        Parameters
        ----------
        index :

        Returns
        -------
        target_index

        """
        return self.src_to_trg.get(index, self.unknown_ind)


class RegionIndexMappingLobe(object):
    """
    Class wrapper for a region index mapping to lobes.

    This allows mapping of each index of source file to an index in lobe file.

    """

    def __init__(
        self,
        orig_annot_file: os.PathLike,
        lobe_annot_file: os.PathLike,
        hemisphere: Hemisphere,
    ):
        # read in the region mapping for each index in annotation file
        region_mapping, _, region_names = nibabel.freesurfer.io.read_annot(
            orig_annot_file
        )
        region_mapping[region_mapping == -1] = 0  # Unknown regions in hemispheres

        # read in the lobe names corresponding to each index in the annotation file
        lobe_mapping, _, lobe_names = nibabel.freesurfer.io.read_annot(lobe_annot_file)
        lobe_mapping[lobe_mapping == -1] = 0

        region_to_lobe_dict = {}

        # get an index for each unique region name
        reginds = []
        for idx, regname in np.unique(region_names):
            regind = np.where(region_names == regname)[0]
            reginds.append(regind)

        # get corresponding lobes names for those indices
        for regind in reginds:
            # get corresponding index in lobe names for this
            lobeind = lobe_mapping[regind]
            lobename = lobe_names[lobeind]

            region_to_lobe_dict[region_names[regind]] = lobename

        self.src_to_trg = region_to_lobe_dict

    def source_to_target(self, index):
        """
        Convert source index on the Source LUT to the Target LUT index.

        Parameters
        ----------
        index :

        Returns
        -------
        target_index

        """
        return self.src_to_trg.get(index, self.unknown_ind)


class ColorLut(object):
    """
    Class wrapper for the color lookup table.

    Each column represents:
    id, name, R, G, B, A, shortname

    """

    def __init__(self, filename: os.PathLike):
        table = np.genfromtxt(os.fspath(filename), encoding="latin1", dtype=None)
        if len(table.dtype) == 6:
            # id name R G B A
            self.inds = table[table.dtype.names[0]]
            self.names = table[table.dtype.names[1]].astype("U")
            self.r = table[table.dtype.names[2]]
            self.g = table[table.dtype.names[3]]
            self.b = table[table.dtype.names[4]]
            self.a = table[table.dtype.names[5]]
            self.shortnames = np.zeros(len(self.names), dtype="U")

        elif len(table.dtype) == 7:
            # id shortname name R G B A
            self.inds = table[table.dtype.names[0]]
            self.shortnames = table[table.dtype.names[1]].astype("U")
            self.names = table[table.dtype.names[2]].astype("U")
            self.r = table[table.dtype.names[3]]
            self.g = table[table.dtype.names[4]]
            self.b = table[table.dtype.names[5]]
            self.a = table[table.dtype.names[6]]
