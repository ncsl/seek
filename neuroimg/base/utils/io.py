import numpy as np
import scipy.io
from mne_bids.utils import _to_tsv
from collections import OrderedDict
from neuroimg.base.utils.utils import MatReader


def save_organized_elecdict_astsv(elecdict, output_fpath, size=None):
    """Save organized electrode dict coordinates as a tsv file."""
    x, y, z, names = list(), list(), list(), list()
    for ch, ch_coord in elecdict.items():
        x.append(ch_coord[0])
        y.append(ch_coord[1])
        z.append(ch_coord[2])
        names.append(ch)
    if size is None:
        sizes = ["n/a"] * len(names)
    else:
        sizes = [size] * len(names)

    data = OrderedDict(
        [("name", names), ("x", x), ("y", y), ("z", z), ("size", sizes),]
    )
    _to_tsv(data, output_fpath)


def save_organized_elecdict_asmat(elecdict, outputfilepath):
    """
    Save centroids as .mat file with attributes eleclabels, which stores
    channel name, electrode type, and depth/grid/strip/other and elecmatrix,
    which stores the centroid coordinates.

    Parameters
    ----------
        elecdict: dict(str: dict(str: np.array))
            Dictionary of channels grouped by electrode.

        outputfilepath: str
            Filepath to save .mat file with annotated electrode labels.
    """
    eleclabels = []
    elecmatrix = []
    # for elec in elecdict.keys():
    for ch_name, coord in elecdict.items():
        label = [[ch_name.strip()], "stereo", "depth"]
        eleclabels.append(label)
        elecmatrix.append(coord)
    mat = {"eleclabels": eleclabels, "elecmatrix": elecmatrix}
    scipy.io.savemat(outputfilepath, mat)


def load_elecs_data(elecfile):
    """
    Load each brain image scan as a NiBabel image object.

    Parameters
    ----------
        elecfile: str
            Path to space-delimited text file of contact labels and contact
            coordinates in mm space.

    Returns
    -------
        eleccoord_mm: dict(str: ndarray)
            Dictionary of contact coordinates in mm space. Maps contact labels
            to contact coordinates, stored as 1x3 numpy arrays.
    """

    eleccoords_mm = {}

    if elecfile.endswith(".txt"):
        with open(elecfile) as f:
            for l in f:
                row = l.split()
                if len(row) == 4:
                    eleccoords_mm[row[0]] = np.array(list(map(float, row[1:])))
                elif len(row) == 6:
                    eleccoords_mm[row[1]] = np.array(list(map(float, row[2:5])))
                else:
                    raise ValueError("Unrecognized electrode coordinate text format")
    else:
        matreader = MatReader()
        data = matreader.loadmat(elecfile)

        eleclabels = data["eleclabels"]
        elecmatrix = data["elecmatrix"]
        print(f"Electrode matrix shape: {elecmatrix.shape}")

        for i in range(len(eleclabels)):
            eleccoords_mm[eleclabels[i][0].strip()] = elecmatrix[i]

    print(f"Electrode labels: {eleccoords_mm.keys()}")

    return eleccoords_mm
