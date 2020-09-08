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

.. image:: https://codeclimate.com/github/ncsl/seek/maintainability
   :target: https://api.codeclimate.com/v1/badges/2c7d5910e89350b967c8/maintainability
   :alt: Code Climate

.. image:: https://img.shields.io/github/repo-size/ncsl/seek
   :target: https://img.shields.io/github/repo-size/ncsl/seek
   :alt: GitHub repo size

.. image:: https://zenodo.org/badge/160566959.svg
   :target: https://zenodo.org/badge/latestdoi/160566959
   :alt: DOI

.. image:: https://api.netlify.com/api/v1/badges/d36d01d2-319a-4e0d-b84f-1d5b4133d5f8/deploy-status
   :target: https://app.netlify.com/sites/elated-almeida-a25d64/deploys
   :alt: Netlify Status

.. image:: https://badges.gitter.im/ncsl/seek.svg
   :target: https://gitter.im/ncsl/seek?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
   :alt: Gitter


This repo describes efforts to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT,
and iEEG data (ECoG, or SEEG). At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry,
regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.
Localizing SEEG electrodes requires at least two contacts on each electrode to initialize the algorithm.
These can be say the deepest 2 contacts, or the entry point and target point (e.g. first and last contact on the electrode).

For ECoG data, we do not explicitly have a process outlined, but these are significantly easier since grids can
be easily interpolated. See `Fieltrip Toolbox`_.

For incorporation of DTI data, see `ndmeg <https://github.com/neurodata/ndmg>`_.

Documentation
-------------

* Link to documentation <>
* For a detailed description of the SEEK workflow of contact localization, specifically localizing the 2 points per electrode, see :doc:`localization guide <./localization_guide>`
* For a detailed description of the overall SEEK workflow, see :ref:`PIPELINEDESCRIPTION`.


Setup and Installation
----------------------

See :doc:`INSTALLATION GUIDE <./INSTALLATION.rst>`_ for full instructions. SEEK uses the Snakemake_
workflow management system to create the different workflows. We chose this because
it is easy to run individual workflows, as well as an entire workflow from the command line.
The full repository is set up similar to the `cookiecutter` Snakemake file: `cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow`.

The recommended installation is via Docker_. See here for instructions on running workflows in the container are shown here below:

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

Creating persistent volumes in Docker
-------------------------------------

If one wants to make a persistent data volume that reflects changes in the Docker container running Snakemake workflows, 
then one can just make a ``data/`` directory inside this repository. Then add in sourcedata. This
directory serves as the BIDS root of the workflows.


Data Organization
-----------------

We use BIDS. 
See https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy

Before data is converted to BIDS in ``seek/pipeline/01-prep`` pipeline, 
then ``sourcedata/`` should contain a semi-structured format of the neuroimaging data that will
be put through the workflow.

**sourcedata/**

.. code-block::

   /{subject}/
       - premri/*.dcm
       - posmri/*.dcm
       - postct/*.dcm


Development
===========

Seek was created and is maintained by `Adam Li <https://adam2392.github.io>`_. It is also maintained and contributed by
Christopher Coogan and other researchers in the NCSL and Crone lab. Contributions are more than welcome so feel free to contact me, open an issue or submit a pull request! See the
:doc:`contribution guide <./doc/contributing>`.

To see the code or report a bug, please visit the `GitHub repository <https://github.com/ncsl/seek>`_.

Note that this program is provided with NO WARRANTY OF ANY KIND. If you can, always double check the results with a human researcher, or clinician.

Pipeline Process Visualized
============================

`DAG of Pipeline in Snakemake <seek/neuroimg/pipeline/dag_neuroimaging_pipeline_reconstruction.pdf>`_

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