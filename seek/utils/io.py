import bz2
import contextlib
import json
import pickle

import numpy as np
import scipy.io


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
        # print(f"Electrode matrix shape: {elecmatrix.shape}")

        for i in range(len(eleclabels)):
            eleccoords_mm[eleclabels[i][0].strip()] = elecmatrix[i]

    # print(f"Electrode labels: {eleccoords_mm.keys()}")

    return eleccoords_mm


class MatReader:
    """
    Object to read mat files into a nested dictionary if need be.
    Helps keep structure from matlab similar to what is used in python.
    """

    def __init__(self, filename=None):
        self.filename = filename

    def loadmat(self, filename):
        """
        this function should be called instead of direct spio.loadmat
        as it cures the problem of not properly recovering python dictionaries
        from mat files. It calls the function check keys to cure all entries
        which are still mat-objects
        """
        data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
        return self._check_keys(data)

    def _check_keys(self, dict):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in dict:
            if isinstance(dict[key], scipy.io.matlab.mio5_params.mat_struct):
                dict[key] = self._todict(dict[key])
        return dict

    def _todict(self, matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        dict = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                dict[strg] = self._todict(elem)
            elif isinstance(elem, np.ndarray):
                dict[strg] = self._tolist(elem)
            else:
                dict[strg] = elem
        return dict

    def _tolist(self, ndarray):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, scipy.io.matlab.mio5_params.mat_struct):
                elem_list.append(self._todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(self._tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    def convertMatToJSON(self, matData, fileName):
        jsonData = {}

        for key in matData.keys():
            if (type(matData[key])) is np.ndarray:
                serializedData = pickle.dumps(
                    matData[key], protocol=0
                )  # protocol 0 is printable ASCII
                jsonData[key] = serializedData
            else:
                jsonData[key] = matData[key]

        with contextlib.closing(bz2.BZ2File(fileName, "wb")) as f:
            json.dump(jsonData, f)
