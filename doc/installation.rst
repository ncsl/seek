.. _installation:

INSTALLATION GUIDE
==================

The best way to install is via a Docker installation.

Docker Installation
-------------------

To run the SEEK pipeline in Docker, first follow instructions to install `Docker <https://docs.docker.com/get-docker/>`_.

**NOTE: You will need approximately at least 8-9 GB free disc space to run the Docker container.**

To setup the container in your system:

.. code-block::

   # build the composition in `docker-compose.yml`
   docker-compose up --build

   # run the container
   docker-compose up


Now if you type in ``docker container ls``\ , you should see the corresponding container.

.. code-block::

   # turn recipe to image
   docker build <image_container_name>

   # turn image to containeer
   docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash

Manual Installation (Not Recommended; See Docker)
-------------------------------------------------

For purposes of documentation and transparency to users, we outline here the manual installation process SEEK can take.
To install the SEEK pipeline manually, one must install the necessary python runtimes, as well as the necessary 3rd party
softwares. 

Python Installations
^^^^^^^^^^^^^^^^^^^^

There are a couple of tools that you need to install in your system before everything is working. You ar recommended to use a Linux based OS. 
Follow links and tutorials on each respective tool to install. Preferably this is done via Docker, or Singularity, but if not, then:

Anaconda and Python3.6+ :


::

    Conda (https://docs.anaconda.com/anaconda/install/)

      This is mainly necessary to run snakemake, and any Python wrapper code.

    .. code-block::

        conda env create -f environment.yml --name=seek
        source activate seek
        conda install sphinx sphinx-gallery sphinx_bootstrap_theme numpydoc black pytest pytest-cov coverage codespell pydocstyle
        pip install coverage-badge anybadge
        # dev versions of mne-python, mne-bids
        pip install --upgrade --no-deps https://api.github.com/repos/mne-tools/mne-python/zipball/master
        pip install --upgrade https://api.github.com/repos/mne-tools/mne-bids/zipball/master

::

    Pip and setup.py install

    .. code-block::

        # run installation via setup.py
        make inplace-all

        # install testing functionality
        make install-tests

        make check
        make test


Pipeline Installations (3rd Party Modules to Install)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Octave

Runs open-source MATLAB-like functions. This runs various scripts for converting output files to object files for rendering visualizations.
Follow: https://www.gnu.org/software/octave/#install

.. code-block::

   brew install octave


#. Gawk_

Runs command line tools.

#. Blender_

Allows nice 3D mesh creations

#. Reconstruction (Freesurfer_)
This step is necessary to generate a parcellation and surface reconstruction of the patient's brain.
The general requirements is just a Linux, or OSX computer with enough RAM.
Currently, this repo is designed to work with FreeSurfer.

#. Coregistration (`FSL Flirt`_)

This step is necessary to map different imaging sessions together. Specifically, for this pipeline, we need it to map CT images to T1 MRI
Note that as of 2019, installation still requires Python2, which should come in any Linux distribution.

     .. code-block::

          python2 <run_installer>

#. Utility (MRTrix3_)

#. SPM_ (preferably 12):

#. Contact-Localization Software (FieldTripToolbox, Img_Pipe, MATLAB)

   * FieldTripToolbox_

#. `ACPC Auto Detection (V2) <https://www.nitrc.org/projects/art/>`:


.. _Gawk: https://brewinstall.org/Install-gawk-on-Mac-with-Brew/
.. _Blender: https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/
.. _Freesurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
.. _FSL Flirt: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/
.. _MRTrix3: https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html
.. _SPM: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
.. _FieldTripToolbox: http://www.fieldtriptoolbox.org/download/