import os
import sys
from distutils.core import setup

import numpy
from setuptools import find_packages

"""
To re-setup: 

    python setup.py sdist bdist_wheel

    pip install -r requirements.txt --process-dependency-links

To test on test pypi:
    
    twine upload --repository testpypi dist/*
    
    # test upload
    pip install -i https://test.pypi.org/simple/ --no-deps neuroimgpipe

    twine upload dist/* 
"""

PACKAGE_NAME = "neuroimgpipe"
with open(os.path.join('seek', '__init__.py'), 'r') as fid:
    for line in (line.strip() for line in fid):
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip('\'').strip('"')
            break
if version is None:
    raise RuntimeError('Could not determine version')
DESCRIPTION = "Neuroimaging Pipeline software for easily generating anatomical interpretations of iEEG data."
URL = "https://github.com/adam2392/neuroimg_pipeline/"
MINIMUM_PYTHON_VERSION = 3, 6  # Minimum of Python 3.6
REQUIRED_PACKAGES = [
    "numpy>=1.14.5",
    "scipy>=1.1.0",
    "pandas>=1.0.0",
    "pybids>=0.10",
    "pybv>=0.2.0",
    "joblib>=0.14",
    "natsort",
    "tqdm",
    "xlrd",
    "nibabel",
    "dicom2nifti",
    "snakemake",
    "matplotlib>=3.2.1",
    "seaborn",
    "mne>=0.20.0",
    "mne-bids>=0.4",
]
CLASSIFICATION_OF_PACKAGE = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    "Development Status :: 3 - Alpha",
    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: Implementation",
    "Natural Language :: English",
]
AUTHORS = [
    "Adam Li",
    "Chester Huynh",
    "Christopher Coogan",
]


def check_python_version():
    """Exit when the Python version is too low."""
    if sys.version_info < MINIMUM_PYTHON_VERSION:
        sys.exit("Python {}.{}+ is required.".format(*MINIMUM_PYTHON_VERSION))


check_python_version()

setup(
    name=PACKAGE_NAME,
    version=version,
    description=DESCRIPTION,
    author=AUTHORS,
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url=URL,
    license="GNU General Public License (GPL)",
    packages=find_packages(exclude=["tests"]),
    project_urls={
        "Documentation": "https://github.com/adam2392/neuroimg_pipeline/docs/",
        "Source": URL,
        "Tracker": "https://github.com/adam2392/neuroimg_pipeline/issues",
    },
    include_dirs=[numpy.get_include()],
    install_requires=REQUIRED_PACKAGES,
    include_package_data=True,
    classifiers=CLASSIFICATION_OF_PACKAGE,
)
