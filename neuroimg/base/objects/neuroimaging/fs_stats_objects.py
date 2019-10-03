import numpy as np


class FSCorticalStats:
    """

    # StructName NumVert SurfArea GrayVol ThickAvg ThickStd MeanCurv GausCurv FoldInd CurvInd
    """

    def __init__(self, filepath):
        self.table = np.genfromtxt(filepath, encoding="latin1", dtype=None)

    @property
    def numverts(self):
        return self.table[self.table.dtype.names[1]]

    @property
    def names(self):
        return self.table[self.table.dtype.names[0]].astype("U")

    @property
    def surfareas(self):
        return self.table[self.table.dtype.names[2]]

    @property
    def volumes(self):
        return self.table[self.table.dtype.names[3]]


class FSSegmentationStats:
    """
     Index SegId NVoxels Volume_mm3 StructName normMean normStdDev normMin normMax normRange
    """

    def __init__(self, filepath):
        self.table = np.genfromtxt(filepath, encoding="latin1", dtype=None)

    @property
    def indices(self):
        # begin indices at 0
        return self.table[self.table.dtype.names[0]] - 1

    @property
    def numvoxels(self):
        return self.table[self.table.dtype.names[1]]

    @property
    def volumes(self):
        return self.table[self.table.dtype.names[2]]

    @property
    def names(self):
        return self.table[self.table.dtype.names[3]].astype("U")
