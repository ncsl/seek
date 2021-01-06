=======================================================
SEEK Pipeline (Stereotactic ElectroEncephalography Kit)
=======================================================

.. image:: https://circleci.com/gh/ncsl/seek.svg?style=svg
   :target: https://circleci.com/gh/ncsl/seek
   :alt: CircleCI

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

.. image:: https://badges.gitter.im/ncsl/seek.svg
   :target: https://gitter.im/ncsl/seek?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
   :alt: Gitter

This repo describes efforts to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT,
and iEEG data (ECoG, or SEEG). For ECoG data, we do not explicitly have a process outlined, but these are significantly easier since grids can
be easily interpolated. See `Fieldtrip Toolbox`_.

We build and keep up-to-date a number of Dockerfiles to run various workflows:

- acpcdetect:
.. image:: https://img.shields.io/docker/image-size/neuroseek/acpcdetect
    :target: https://hub.docker.com/repository/docker/neuroseek/acpcdetect
    :alt: Docker Image Size (tag)

- blender:
.. image:: https://img.shields.io/docker/image-size/neuroseek/blender
    :target: https://hub.docker.com/repository/docker/neuroseek/blender
    :alt: Docker Image Size (tag)

- FreeSurfer7+ with MRTrix3:
.. image:: https://img.shields.io/docker/image-size/neuroseek/freesurfer7-with-mrtrix3
    :target: https://hub.docker.com/repository/docker/neuroseek/freesurfer7-with-mrtrix3
    :alt: Docker Image Size (tag)


Note: For FSL, we use a 3rd party docker image.

Documentation
-------------

* For a detailed description of the overall SEEK workflow, see `workflow documentation <https://github.com/ncsl/seek/blob/master/workflow/documentation.md>`_.
* For a detailed description of the SEEK workflow of contact localization, specifically localizing the 2 points per electrode, see :doc:`localization guide <./tutorials/localization_guide>`_

For a description of the visualization engine, see: https://github.com/cronelab/ReconstructionVisualizer

Setup and Installation
----------------------

See `INSTALLATION GUIDE <installation.rst>`_ for full instructions. SEEK uses the Snakemake_
workflow management system to create the different workflows. We chose this because
it is easy to run individual workflows, as well as an entire workflow from the command line.
The full repository is set up similar to the `cookiecutter` Snakemake file: `cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow`.

The recommended installation is via Docker_. See here for instructions on running workflows in the container are shown here below:

Data Organization
-----------------

We use BIDS. See https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy

Before data is converted to BIDS in ``seek/pipeline/01-prep`` pipeline,
then ``sourcedata/`` should contain a semi-structured format of the neuroimaging data that will
be put through the workflow.

**sourcedata/**

.. code-block::

   /{subject}/
       - premri/*.dcm
       - posmri/*.dcm
       - postct/*.dcm


Docker Installation
-------------------

The Docker containers sit on `Docker Hub`_, specifically at `https://hub.docker.com/r/neuroseek/seek <https://hub.docker.com/r/neuroseek/seek>`_.

Setup: Note that the docker container names are:

.. code-block::

   - seek_reconstruction
   - seek_localization  # tbd
   - seek_visualization  # tbd


To setup the container in your system:

.. code-block::

   docker-compose up --build

::

    Running workflows within the container:
    In another terminal, one can run the pipeline commands in terminal.

    .. code-block::

       # run container and mount data directories
       docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash

::

    Running workflows **using** the container:

    .. code-block::

       # run snakemake using the containers
       snakemake <rule_name> --use-singularity

For running individual pipelines, see `INSTALLATION GUIDE <INSTALLATION.md>`_.

Running workflows using Docker and Snakemake
--------------------------------------------
To run snakemake workflows using Docker, we have implemented various ``Makefile`` recipes.

See:

* snakemake-all
* recon
* prep-localization
* coregistration
* prep-viz

``You will need to alter the bind paths to your specific BIDS root directory and SEEK repository directory``.

Creating persistent volumes in Docker
-------------------------------------

If one wants to make a persistent data volume that reflects changes in the Docker container running Snakemake workflows, 
then one can just make a ``data/`` directory inside this repository. Then add in sourcedata. This
directory serves as the BIDS root of the workflows.


Development
===========

Seek was created and is maintained by `Adam Li <https://adam2392.github.io>`_. It is also maintained and contributed by
Christopher Coogan and other researchers in the NCSL and Crone lab. Contributions are more than welcome so feel free to contact me, open an issue or submit a pull request! See the
:doc:`contribution guide <./doc/contributing>`.

To report a bug, please visit the `GitHub repository <https://github.com/ncsl/seek/issues/>`_.

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


.. _Docker: https://www.docker.com/
.. _Docker Hub: https://hub.docker.com/
.. _FieldTrip Toolbox: http://www.fieldtriptoolbox.org/tutorial/human_ecog/
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/

FAQ
===
1. For incorporation of DTI data, see `ndmeg <https://github.com/neurodata/ndmg>`_.