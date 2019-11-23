Install
=======


Below we assume you have the default Python environment already configured on
your computer and you intend to install ``neuroimgpipe`` inside of it.  If you want
to create and work with Python virtual environments, please follow instructions
on `venv <https://docs.python.org/3/library/venv.html>`_ and `virtual
environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.

First, make sure you have the latest version of ``pip`` (the Python package manager)
installed. If you do not, refer to the `Pip documentation
<https://pip.pypa.io/en/stable/installing/>`_ and install ``pip`` first.

Install the released version
----------------------------

This pipeline is meant for you to adjust data paths and run locally.
You can manually download ``neuroimgpipe`` from
`GitHub <https://github.com//neuroimg_pipeline/releases>`_


Python package dependencies
---------------------------
neuroimgpipe requires the following packages:

- numpy
- pandas
- scipy
- snakemake
- skimage


Hardware requirements
---------------------
`neuroimgpipe` package requires only a standard computer with enough RAM to support the in-memory operations.

OS Requirements
---------------
This package is supported for *Linux* and *macOS*. However, the package has been tested on the following systems:

- Linux: N/A
- macOS: N/A