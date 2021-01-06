"""Neuroimaging Pipeline software for easily generating anatomical interpretations of iEEG data."""

__name__ = "seek"
__version__ = "0.1.0"

from .format import bids_validate, convert_img_to_bids
from .utils import BidsRoot, FREESURFER_HOME
