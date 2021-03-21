=======================================================
SEEK-Pipeline (Stereotactic ElectroEncephalography Kit)
=======================================================

.. image:: https://circleci.com/gh/ncsl/seek.svg?style=svg
   :target: https://circleci.com/gh/ncsl/seek
   :alt: CircleCI

.. image:: https://github.com/ncsl/seek/workflows/.github/workflows/main.yml/badge.svg
    :target: https://github.com/ncsl/seek/actions/
    :alt: GitHub Actions

.. image:: https://travis-ci.com/ncsl/seek.svg?token=6sshyCajdyLy6EhT8YAq&branch=master
   :target: https://travis-ci.com/ncsl/seek
   :alt: Build Status

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Code style: black

.. image:: https://img.shields.io/github/repo-size/ncsl/seek
   :target: https://img.shields.io/github/repo-size/ncsl/seek
   :alt: GitHub repo size

.. image:: https://zenodo.org/badge/160566959.svg
   :target: https://zenodo.org/badge/latestdoi/160566959
   :alt: DOI

.. image:: https://img.shields.io/badge/snakemake-≥5.27.4-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake

Stereotactic EEG Kit Pipeline (SEEK-Pipeline) is a `Snakemake`_ pipeline workflow for running T1w reconstruction (via `FreeSurfer`_),
CT -> T1w coregistration, BIDS-formatting and mesh generation (via `Blender`_).

A few notes
-----------

- Although we primarily focus on sEEG implantations, this pipeline can hypothetically work well with ECoG data.
- Although the pipeline is mostly automated via a `config file <https://github.com/ncsl/seek/blob/master/config/localconfig.yaml>`_, the actual contact localization of the iEEG electrodes on a CT/T1w image requires a GUI. We currently use `Fieldtrip Toolbox`_.

Pipeline containerization
-------------------------

In order to create a 100% reproducible workflow that also limits the need to install various complex software
on your systems, we utilize Docker and Singularity that just "works" on snakemake.
We build and keep up-to-date a number of Dockerfiles to run various workflows:

- acpcdetect: |acpcdetect|
- blender: |blender|
- FreeSurfer7+ with MRTrix3: |freesurfer|

.. |acpcdetect| image:: https://img.shields.io/docker/image-size/neuroseek/acpcdetect
    :target: https://hub.docker.com/repository/docker/neuroseek/acpcdetect
    :alt: Docker Image Size (tag)
.. |blender| image:: https://img.shields.io/docker/image-size/neuroseek/blender
    :target: https://hub.docker.com/repository/docker/neuroseek/blender
    :alt: Docker Image Size (tag)
.. |freesurfer| image:: https://img.shields.io/docker/image-size/neuroseek/freesurfer7-with-mrtrix3
    :target: https://hub.docker.com/repository/docker/neuroseek/freesurfer7-with-mrtrix3
    :alt: Docker Image Size (tag)

Note: For FSL, we use a 3rd party docker image.

Documentation
-------------
SEEK is comprised of more then just a pipeline. We also have more explicit documentation for exactly
how to localize electrodes using the localization GUI we recommend. We also have a 3D web-based visualization
engine that can render dynamic 3/4-D visualizations of the sEEG data.

To see the entire documentation, see http://neuroseek.azurewebsites.net/docs/seek/

* For a detailed description of the overall SEEK workflow, see `workflow documentation <https://github.com/ncsl/seek/blob/master/workflow/documentation.md>`_.
* For a detailed description of the SEEK workflow of contact localization, specifically localizing the 2 points per electrode, see `localization guide <http://neuroseek.azurewebsites.net/docs/localize/>`_
* For the visualization engine, see: https://github.com/cronelab/ReconstructionVisualizer

For detailed setup, installation and usage instructions, see documentation.

Installation
------------

SEEK uses the Snakemake_ workflow management system to create different workflows. We chose this because
it is easy to run individual workflows, as well as an entire workflow from the command line. To install the 
pipeline from github:

.. code-block:: bash

    # clone repository locally
    git clone https://github.com/ncsl/seek
    python3.8 -m venv .venv
    pipenv install

You will also need Docker_ and Singularity_ and some configuration setup to be able to successfully run the workflows. 
See documentation for full instructions. 

Running workflows using Docker and Snakemake
--------------------------------------------
To run snakemake workflows using Docker, we have implemented various ``Makefile`` recipes.
First, you need to change directory to the correct workflow and then run snakemake with some
arguments to bind directories to the singularity container. For example

.. code-block:: bash

    snakemake --cores 1 --use-singularity --singularity-args "--bind ~/hdd/epilepsy_bids/,~/Documents/seek/";

where, you can alter the cores used, and also bind various directories.
See ``Makefile`` and documentation for more details on the following pipelines:

* snakemake-all
* recon
* prep-localization
* coregistration
* prep-viz

``You will need to alter the bind paths to your specific BIDS root directory and SEEK repository directory``.

Development
===========

Seek was created and is maintained by `Adam Li <https://adam2392.github.io>`_. It is also maintained and contributed by
`Christopher Coogan <https://github.com/TheBrainChain>`_ and other researchers in the NCSL and Crone lab. Contributions are more than welcome so feel free to contact me, open an issue or submit a pull request! See the
:doc:`contribution guide <./doc/contributing>`. To report a bug, please visit the `GitHub repository <https://github.com/ncsl/seek/issues/>`_.

Note that this program is provided with NO WARRANTY OF ANY KIND. If you can, always double check the results with a human researcher, or clinician.

How to cite SEEK?
=================

If you want to cite Seek, please use the Zenodo for the repository.

Acknowledgement
===============

Several functions of Seek essentially make use of existing software packages for neuroimaging analysis, including:

- `Recon-all (FreeSurfer) <https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all>`_
- `FSL Flirt <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_
- `MRTrix3 <http://www.mrtrix.org/>`_
- `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_

.. _FreeSurfer: https://surfer.nmr.mgh.harvard.edu/
.. _Blender: https://www.blender.org/
.. _Docker: https://www.docker.com/
.. _Docker Hub: https://hub.docker.com/
.. _FieldTrip Toolbox: http://www.fieldtriptoolbox.org/tutorial/human_ecog/
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/

FAQ
===
1. For incorporation of DTI data, see `ndmeg <https://github.com/neurodata/ndmg>`_.
