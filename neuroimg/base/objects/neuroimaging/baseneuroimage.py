import enum
import os
import re
import subprocess
import tempfile
from typing import List, Optional

import nibabel
import numpy as np

"""
Stores the base classes for neuroimaging interaction.

- Hemisphere: an enum object to denote which side of the brain 
- ColorLut: wrapper to access the color LUT.
- RegionIndexMapping: wrapper to map region indices from one LUT -> another LUT

"""

SUBCORTICAL_REG_INDS = [8, 10, 11, 12, 13, 16, 17, 18, 26, 47, 49, 50, 51, 52, 53, 54, 58]
FS_LUT_LH_SHIFT = 1000
FS_LUT_RH_SHIFT = 2000

class Hemisphere(enum.Enum):
    rh = 'rh'
    lh = 'lh'

class RegionIndexMapping:
    """
    Class wrapper for a region index mapping.

    This maps each index of the source file to an index in the corresponding look up table.
    """

    def __init__(self, color_lut_src_file: os.PathLike, color_lut_trg_file: os.PathLike):
        """
        :param color_lut_src_file: The source lookup table
        :param color_lut_trg_file: The target's lookup table.
        """
        self.src_table = ColorLut(color_lut_src_file)
        self.trg_table = ColorLut(color_lut_trg_file)

        names_to_trg = dict(zip(self.trg_table.names, self.trg_table.inds))

        self.src_to_trg = dict()
        for src_ind, src_name in zip(self.src_table.inds, self.src_table.names):
            trg_ind = names_to_trg.get(src_name, None)
            if trg_ind is not None:
                self.src_to_trg[src_ind] = trg_ind

        self.unknown_ind = names_to_trg.get('Unknown', 0)  # zero as the default unknown area

    def source_to_target(self, index):
        return self.src_to_trg.get(index, self.unknown_ind)


def pial_to_verts_and_triangs(pial_surf) -> (np.ndarray, np.ndarray):
    """
    Function to convert pial surface file to vertices and triangles

    :param pial_surf:
    :return:
    """
    tmpdir = tempfile.TemporaryDirectory()
    pial_asc = os.path.join(tmpdir.name, os.path.basename(pial_surf + ".asc"))
    subprocess.run(['mris_convert', pial_surf, pial_asc])

    with open(pial_asc, 'r') as f:
        f.readline()
        nverts, ntriangs = [int(n) for n in f.readline().strip().split(' ')]

    vertices = np.genfromtxt(pial_asc, dtype=float, skip_header=2, skip_footer=ntriangs, usecols=(0, 1, 2))
    triangles = np.genfromtxt(pial_asc, dtype=int, skip_header=2 + nverts, usecols=(0, 1, 2))
    assert vertices.shape == (nverts, 3)
    assert triangles.shape == (ntriangs, 3)

    completed_process = subprocess.run(["mris_info", pial_surf], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    mris_info = completed_process.stdout.decode('ascii')
    c_ras_list = extract_vector(mris_info, "c_(ras)")
    assert c_ras_list is not None
    vertices[:, 0:3] += np.array(c_ras_list)

    return vertices, triangles


def extract_vector(string: str, name: str) -> Optional[List[float]]:
    r"""
    Extract numerical vector from a block of text. The vector has to be on a single line with the format:
    <name> : (x0, x1, x2 [,...])
    If the vector in the correct format is missing, return None.

    >>> extract_vector("EXAMPLE\na: (1.0, 2.0, 3.0)\nb: (0.0, 0.0, 0.0)", "a")
    [1.0, 2.0, 3.0]

    >>> extract_vector("EMPTY", "a") is None
    True
    """

    for line in string.split("\n"):
        match = re.match(r"""^\s*
                         (.+?)              # name
                         \s*:\s*            # separator
                         \(([0-9.,\s-]+)\)   # vector: (x0, x1, ....)
                         \s*$""",
                         line, re.X)

        if match and match.group(1) == name:
            try:
                vector = [float(x) for x in match.group(2).split(",")]
                return vector
            except ValueError:
                pass

    return None


def read_cortical_region_mapping(label_direc: os.PathLike, hemisphere: Hemisphere, fs_to_conn: RegionIndexMapping) \
        -> np.ndarray:
    """
    Reads the cortical region mapping file.

    :param label_direc: Where the annotation label directory is
    :param hemisphere: (Hemisphere) enumerator
    :param fs_to_conn: (RegionIndexMapping)
    :return:
    """
    filename = os.path.join(label_direc, hemisphere.value + ".aparc.annot")
    region_mapping, _, _ = nibabel.freesurfer.io.read_annot(filename)

    region_mapping[region_mapping == -1] = 0  # Unknown regions in hemispheres

    # $FREESURFER_HOME/FreeSurferColorLUT.txt describes the shift
    if hemisphere == Hemisphere.lh:
        region_mapping += FS_LUT_LH_SHIFT
    else:
        region_mapping += FS_LUT_RH_SHIFT

    fs_to_conn_fun = np.vectorize(lambda n: fs_to_conn.source_to_target(n))
    region_mapping = fs_to_conn_fun(region_mapping)

    return region_mapping


class ColorLut:
    """
    Class wrapper for the color lookup table.

    each column represents:
    id, name, R, G, B, A, shortname
    """

    def __init__(self, filename: os.PathLike):
        table = np.genfromtxt(os.fspath(filename), dtype=None)

        if len(table.dtype) == 6:
            # id name R G B A
            self.inds = table[table.dtype.names[0]]
            self.names = table[table.dtype.names[1]].astype('U')
            self.r = table[table.dtype.names[2]]
            self.g = table[table.dtype.names[3]]
            self.b = table[table.dtype.names[4]]
            self.a = table[table.dtype.names[5]]
            self.shortnames = np.zeros(len(self.names), dtype='U')

        elif len(table.dtype) == 7:
            # id shortname name R G B A
            self.inds = table[table.dtype.names[0]]
            self.shortnames = table[table.dtype.names[1]].astype('U')
            self.names = table[table.dtype.names[2]].astype('U')
            self.r = table[table.dtype.names[3]]
            self.g = table[table.dtype.names[4]]
            self.b = table[table.dtype.names[5]]
            self.a = table[table.dtype.names[6]]

