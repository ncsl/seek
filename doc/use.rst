:orphan:

.. _use:

SEEK-Pipeline Usage Guide
=========================
Using seek requires some minimal working knowledge of the command line, snakemake_ and
Docker_. If you have ever manually ran FreeSurfer, coregistered CT images to T1 images,
localized electrode coordinates and then reformatted files to be BIDS-compliant, then
you'll see that seek attempts to abstract away many of those complexities. In addition,
we generate useful Blender_ objects that can be used in our downstream
`visualization platform <https://github.com/cronelab/ReconstructionVisualizer>`_.

If you have not followed installation instructions yet, please visit :doc:`installation page <installation>`
to fully setup for usage.

Configuration Setup
-------------------
Next you will want to modify the snakemake configuration files. For more
detailed information on how snakemake configuration works, see:
https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html

To configure SEEK, we simply need a list of subjects and a series of paths
that are useful. First, modify the config ``yaml`` file.

.. code-block::

    # modify config yaml
    $ cd seek/
    $ vim config/localconfig.yaml

Notice how you will have to input the path to the root of your BIDS dataset,
corresponding to ``bids_root``. In addition, you can specify the ``session``
BIDS-entity, which will describe how certain files are organized. For more info,
see BIDS_.

Next, modify the ``subjects.tsv`` file to specify which subjects you
want to run workflows on.

.. code-block::

    # modify the subjects table
    $ vim config/subjects.tsv

Running workflows
=================
Now that you have setup everything, you are ready to
run the workflows. You can directly call snakemake `commands <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_
from each of the workflow's directories. Otherwise, each workflow has a
`Makefile` recipe command that can be called.

The workflows can be ran in sequence as described below.
For more detailed description of the pipelines themselves, that are automated, see :doc:`pipeline description <pipeline_description>`

Multiple cores can be used to basically run multiple rules at once. Snakemake will handle 
underneath the hood the rule dependencies.

Setting Environment
-------------------
SEEK requires ``snakemake`` to run. Usually this should be setup in a virtual
environment from :doc:`installation instructions <installation>`. To then
setup environment variables ``SEEKHOME``, which is the path to your seek directory,
one should run

.. code-block::

    # activate your virtual environment (here we assume pipenv)
    pipenv shell
    export SEEKHOME=<path_to_seek_repo>


.. _dicom_to_bids

1. Converting Dicoms to Nifti as BIDS
-------------------------------------

If you already have Nifti files BIDS-formatted, then you can skip to :ref:`t1w_process`.

The raw form of imaging data comes in the form of Dicoms, which must be organized according 
to the :doc:`data setup instructions in installation <installation>`. This step 
would use a containerized and conda environment to run ``heudiconv`` underneath the hood 
to convert ``*.dcm`` files into ``*.nii.gz`` files that are BIDS-compliant. Initially your data 
might look something like this:

.. image:: /_static/dicom_start.png
    :width: 400
    :alt: Starting sourcedata directory

The following command will run this pipeline.

.. code-block::

    make bids

and afterwards your dataset will look something like this:

.. image:: /_static/nifti_bids.png
    :width: 400
    :alt: After BIDS conversion

.. _t1w_process

2. Processing T1w images (FOV Crop, ACPC alignment, Reconstruction)
-------------------------------------------------------------------
This will abstract away FreeSurfer reconstruction commands and various other
commands.

This will organize FreeSurfer output into the ``<bids_root>/derivatives/freesurfer/<subject_id>``
folders and also generate other derivatives.

.. code-block::

    make recon

Prep Localization (prep_localization)
-------------------------------
This will run very simple rules to setup the input for electrode localization.

.. code-block::

    make prep-localize

Electrodes Localization
-----------------------
At this point, SEEK recommends using
``FieldTrip`` toolbox to localize electrodes
on the CT image. We recommend following the
localization tutorial we have built and using
the following matlab script: https://github.com/ncsl/seek/blob/master/workflow/scripts/electrode_localization.m

The output of that script will be an ``*electrodes.tsv`` file.

Label Contacts (contact_labeling)
-------------------------------
This workflow will take the ``*electrodes.tsv`` file generated from localization, and
then generate additional files in different coordinate systems, along with the
``*coordsystem.json`` files, as specified in BIDS_.

.. code-block::

    make label-contacts

Prep Visualization (prep_vizengine)
-------------------------------
This will generate Blender_ objects, which can then be used with SEEK-Viz.

.. code-block::

    make prep-viz

Parallelization (running multiple subjects at once)
===================================================
With the help of snakemake_, every workflow above can be parallelized trivially.
We add a ``cores`` option to each Makefile recipe. For example

.. code-block::

    make recon cores=3

will run the ``recon`` workflow with 3 cores.

.. _snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Docker: https://docs.docker.com/get-docker/
.. _Blender: https://www.blender.org/
.. _BIDS: https://bids-specification.readthedocs.io/