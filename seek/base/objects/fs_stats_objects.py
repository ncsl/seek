import numpy as np


class FSCorticalStats:
    """StructName NumVert SurfArea GrayVol ThickAvg ThickStd MeanCurv GausCurv FoldInd CurvInd."""

    def __init__(self, filepath):
        self.table = np.genfromtxt(filepath, encoding="latin1", dtype=None)

    @property
    def numverts(self):
        """
        Vertices amount for each FreeSurfer Cortical Region.

        Returns
        -------
        numverts :
        """
        return self.table[self.table.dtype.names[1]]

    @property
    def names(self):
        """
        Unicode of the FreeSurfer segmentation names.

        Returns
        -------
        names :

        """
        return self.table[self.table.dtype.names[0]].astype("U")

    @property
    def surfareas(self):
        """
        Surface areas for each FreeSurfer Cortical region.

        Returns
        -------
        surface_areas :
        """
        return self.table[self.table.dtype.names[2]]

    @property
    def volumes(self):
        """
        Volumes of each FreeSurfer segmented region.

        Returns
        -------
        volumes :
        """
        return self.table[self.table.dtype.names[3]]


class FSSegmentationStats:
    """Index SegId NVoxels Volume_mm3 StructName normMean normStdDev normMin normMax normRange."""

    def __init__(self, filepath):
        self.table = np.genfromtxt(filepath, encoding="latin1", dtype=None)

    @property
    def indices(self):
        """
        First column of indices for each FreeSurfer Segmentation statistics.

        Returns
        -------
        indices :
        """
        # begin indices at 0
        return self.table[self.table.dtype.names[0]] - 1

    @property
    def numvoxels(self):
        """
        Voxels for each FreeSurfer segmented region.

        Returns
        -------
        numvoxels :
        """
        return self.table[self.table.dtype.names[1]]

    @property
    def volumes(self):
        """
        Volumes of each FreeSurfer segmented region.

        Returns
        -------
        volumes :
        """
        return self.table[self.table.dtype.names[2]]

    @property
    def names(self):
        """
        Unicode of the FreeSurfer segmentation names.

        Returns
        -------
        names :

        """
        return self.table[self.table.dtype.names[3]].astype("U")
