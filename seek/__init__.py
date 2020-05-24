"""Neuroimaging Pipeline software for easily generating anatomical interpretations of iEEG data."""
from os.path import dirname, basename, isfile
import glob
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

modules = glob.glob(dirname(__file__) + "/*.py")
__all__ = [
    basename(f)[:-3] for f in modules if isfile(f) and not f.endswith("__init__.py")
]
__name__ = "seek"
__version__ = "0.1.0"
