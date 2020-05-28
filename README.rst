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


:raw-html-m2r:`<a href="https://codeclimate.com/github/ncsl/seek/maintainability"><img src="https://api.codeclimate.com/v1/badges/2c7d5910e89350b967c8/maintainability" /></a>`


.. image:: https://img.shields.io/github/repo-size/ncsl/seek
   :target: https://img.shields.io/github/repo-size/ncsl/seek
   :alt: GitHub repo size


.. image:: https://zenodo.org/badge/160566959.svg
   :target: https://zenodo.org/badge/latestdoi/160566959
   :alt: DOI


.. image:: https://images.microbadger.com/badges/version/neuroseek/seek.svg
   :target: https://microbadger.com/images/neuroseek/seek "Get your own version badge on microbadger.com"
   :alt: 


This repo describes Sarma/Crone lab effort to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT,
and iEEG data (ECoG, or SEEG). 

For incorporation of DTI data, see `ndmeg <https://github.com/neurodata/ndmg>`_.


.. raw:: html

   <!-- MarkdownTOC -->


* Features
* Setup and Installation

  * DOCKER

* Creating persistent volumes
  .. code-block::

       - Reconstruction
       - Electrode localization \(Bioimage suite\)

* Data Organization
* Pipeline Description
* Documentation and Testing

  * Pipeline Process Visualized

* References:


.. raw:: html

   <!-- /MarkdownTOC -->



Features
--------


* [ ] Add support for MRICloud running using R-script. Possibly convert to Python script.
* [ ] Create unit and integration tests using pytest that test: pipeline in both snakemake and Python

Setup and Installation
----------------------

See `INSTALLATION GUIDE <INSTALLATION.md>`_. SEEK uses the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ 
workflow management system to create the different workflows. We chose this because
it is easy to run individual workflows, as well as an entire workflow from the command line.

DOCKER
^^^^^^

Setup: Note that the docker container names are:

.. code-block::

   - seek_reconstruction
   - seek_localization  # tbd
   - seek_visualization  # tbd


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


#. 
   Run localization container:

   ..

      docker-compose run localization ./start_bioimagesuite


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



Pipeline Description
--------------------

At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, 
regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

See `PIPELINE GUIDE <PIPELINE_DESCRIPTION.md>`_

Semi-Automated Localizing Electrodes Process

----

Localizing SEEG electrodes requires at least two contacts on each electrode to initialize the algorithm.
These can be say the deepest 2 contacts, or the entry point and target point (e.g. first and last contact on the electrode).

For ECoG data, we do not explicitly have a process outlined, but these are significantly easier since grids can
be easily interpolated.

See `LOCALIZATION_GUIDE <LOCALIZATION_GUIDE.md>`_

Documentation and Testing
-------------------------

See `Testing Guide <TESTING_SETUP.md>`_

Pipeline Process Visualized
^^^^^^^^^^^^^^^^^^^^^^^^^^^

`DAG of Pipeline in Snakemake <seek/neuroimg/pipeline/dag_neuroimaging_pipeline_reconstruction.pdf>`_

References:
-----------

#. Recon-all. FreeSurfer. https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all#References
#. FSL Flirt. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT
#. MRTrix3. http://www.mrtrix.org/
#. Img_pipe. https://github.com/ChangLabUcsf/img_pipe
#. MRICloud. https://mricloud.org/
#. Snakemake. https://snakemake.readthedocs.io/en/stable/
#. FieldTrip Toolbox. http://www.fieldtriptoolbox.org/tutorial/human_ecog/
