import argparse
import os
import sys

import nibabel as nb
import numpy as np
import numpy.linalg as npl
import scipy.io
from nibabel.affines import apply_affine

from neuroimg.base.utils import MatReader
from neuroimg.localize_contacts.electrode_clustering.mask import MaskVolume
from neuroimg.localize_contacts.electrode_clustering.grouping import Cluster, CylindricalGroup
from neuroimg.localize_contacts.electrode_clustering.postprocess import PostProcessor
from neuroimg.localize_contacts.freecog_labeling.utils import convert_fsmesh2mlab, label_elecs

class Test_Localization():
    def test_clustering(self):
        pass