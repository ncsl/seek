SEEK Pipeline (Stereotactic ElectroEncephalography Kit)
-------------------------------------------------------


.. image:: https://circleci.com/gh/ncsl/seek.svg?style=svg
   :target: https://circleci.com/gh/ncsl/seek
   :alt: CircleCI


.. image:: https://travis-ci.com/ncsl/seek.svg?token=6sshyCajdyLy6EhT8YAq&branch=master
   :target: https://travis-ci.com/ncsl/seek
   :alt: Build Status


.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Code style: black


.. image:: https://img.shields.io/github/license/ncsl/seek
   :target: https://img.shields.io/github/license/ncsl/seek
   :alt: GitHub


.. image:: https://img.shields.io/github/last-commit/ncsl/seek
   :target: https://img.shields.io/github/last-commit/ncsl/seek
   :alt: GitHub last commit


.. image:: https://img.shields.io/github/repo-size/ncsl/seek
   :target: https://img.shields.io/github/repo-size/ncsl/seek
   :alt: GitHub repo size


.. image:: https://zenodo.org/badge/160566959.svg
   :target: https://zenodo.org/badge/latestdoi/160566959
   :alt: DOI

.. image:: https://images.microbadger.com/badges/version/neuroseek/seek.svg
   :target: https://microbadger.com/images/neuroseek/seek "Get your own version badge on microbadger.com"
   :alt: version


This repo describes Sarma/Crone lab effort to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT,
and iEEG data (ECoG, or SEEG). 

For incorporation of DTI data, see `ndmeg <https://github.com/neurodata/ndmg>`_.


Features
--------

* [ ] Add support for MRICloud running using R-script. Possibly convert to Python script.
* [ ] Create unit and integration tests using pytest that test: pipeline in both snakemake and Python

Setup and Installation
----------------------

See `installation <./doc/INSTALLATION.rst>`_. SEEK uses the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_
workflow management system to create the different workflows. We chose this because
it is easy to run individual workflows, as well as an entire workflow from the command line.

Docker
------

`DOCKER <https://hub.docker.com/orgs/neuroseek/repositories>`_

Setup: Note that the docker container names are:

.. code-block::

   seek_reconstruction
   seek_localization  # tbd
   seek_visualization  # tbd


To setup the container in your system:

.. code-block::

   docker-compose up --build


In another terminal run the pipeline commands.

.. code-block::

   # turn image to containeer
   docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash


For running individual pipelines, see `INSTALLATION GUIDE <INSTALLATION.md>`_.

Creating persistent volumes
===========================

If one wants to make a persistent data volume that reflects changes in the Docker container running Snakemake workflows, 
then one can just make a ``Data/`` directory inside this repository. Then add in sourcedata. This
directory serves as the BIDS root of the workflows.

Electrode localization (Bioimage suite)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run localization container:

.. code-block::

      docker-compose run localization ./start_bioimagesuite

Data Organization
-----------------

We use BIDS. 
See https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy

Before data is converted to BIDS in ``seek/pipeline/01-prep`` pipeline, 
then ``sourcedata/`` should contain a semi-structured format of the neuroimaging data that will
be put through the workflow.


.. code-block::

    sourcedata/
       /{subject}/
           - premri/*.dcm
           - posmri/*.dcm
           - postct/*.dcm

**Note** Currently this "structured" format of the `sourcedata` is needed in order to run the
Snakemake pipeline, since we have hardcoded where to look for the initial dicoms. In the future,
we are hoping to abstract this away and also allow for multiple scans.

Pipeline Description
--------------------

At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, 
regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

See `PIPELINE GUIDE <PIPELINE_DESCRIPTION.rst>`_

Semi-Automated Localizing Electrodes Process
--------------------------------------------

Localizing SEEG electrodes requires at least two contacts on each electrode to initialize the algorithm.
These can be say the deepest 2 contacts, or the entry point and target point (e.g. first and last contact on the electrode).

For ECoG data, we do not explicitly have a process outlined, but these are significantly easier since grids can
be easily interpolated.

See `LOCALIZATION_GUIDE <LOCALIZATION_GUIDE.rst>`_

Contributing
------------

See `Contribution Guide <contributing.rst>`_. We are always looking for contributors, whether it is an
example or successful use case of the pipeline, extending the pipeline with additional rules, or
contributing documentation.

Pipeline Process Visualized
---------------------------

`DAG of Pipeline in Snakemake <seek/pipeline/dag_neuroimaging_pipeline_reconstruction.pdf>`_

References:
-----------

.. _Gawk: https://brewinstall.org/Install-gawk-on-Mac-with-Brew/
.. _Blender: https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/
.. _Freesurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
.. _FSL Flirt: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/
.. _MRTrix3: https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html
.. _SPM: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
.. _FieldTripToolbox: http://www.fieldtriptoolbox.org/download/