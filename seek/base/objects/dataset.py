import os
import tempfile
from typing import List
from zipfile import ZipFile

import numpy as np


class StructuralDataset:
    """
    Class for a structural type dataset that encompasses brain surfaces.

    Parameters
    ----------
    orientations :
    areas :
    centers :
    cortical :
    names :
    """

    def __init__(
        self,
        orientations: np.ndarray,
        areas: np.ndarray,
        centers: np.ndarray,
        cortical: np.ndarray,
        names: List[str],
    ):
        nregions = len(names)

        # make some assertions on the shape
        assert orientations.shape == (nregions, 3)
        assert areas.shape == (nregions,)
        assert centers.shape == (nregions, 3)
        assert cortical.shape == (nregions,)

        # assign these to class attributes
        self.orientations = orientations
        self.areas = areas
        self.centers = centers
        self.cortical = cortical
        self.names = names

    def save_to_txt_zip(self, filename: os.PathLike):
        """
        Create a zipped file with the necessary .txt files to recreate the entire dataset.

        Includes data for:
        - areas of surface regions
        - average_orientations of each surface region
        - centres in xyz coordinates of each surface
        - index labeling of the cortical

        Parameters
        ----------
        filename :
        """
        # create a temporary directory
        tmpdir = tempfile.TemporaryDirectory()

        file_areas = os.path.join(tmpdir.name, "areas.txt")
        file_orientations = os.path.join(tmpdir.name, "average_orientations.txt")
        file_centres = os.path.join(tmpdir.name, "centres.txt")
        file_cortical = os.path.join(tmpdir.name, "cortical.txt")

        np.savetxt(fname=file_areas, X=self.areas, fmt="%.2f")
        np.savetxt(fname=file_orientations, X=self.orientations, fmt="%.2f %.2f %.2f")
        np.savetxt(fname=file_cortical, X=self.cortical, fmt="%d")

        with open(file_centres, "w") as f:
            for i, name in enumerate(self.names):
                f.write(
                    "%s %.4f %.4f %.4f\n"
                    % (name, self.centers[i, 0], self.centers[i, 1], self.centers[i, 2])
                )

        with ZipFile(filename, "w") as zip_file:
            zip_file.write(file_areas, os.path.basename(file_areas))
            zip_file.write(file_orientations, os.path.basename(file_orientations))
            zip_file.write(file_centres, os.path.basename(file_centres))
            zip_file.write(file_cortical, os.path.basename(file_cortical))
