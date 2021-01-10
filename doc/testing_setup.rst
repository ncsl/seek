:orphan:

.. _testing:

Guide to setup Testing Environment
==================================
Seek is a pipeline manager for the workflows involving SEEG electrode implantations,
T1 MRI and CT scans.

We use `snakemake` module to link the pipelines in a Directed Acyclic Graph (DAG).

Testing Data
------------

We utilize a testing dataset at https://github.com/adam2392/seek-testing-data to
perform unit and integration tests.

Testing Script
--------------
First off, we follow the `black` coding format. In addition, we utilize `pydocstyle`, and `codespell`
checks.

    black seek/*
    black tests/*
    black --check seek/

or

    # the make file has a sequence command to run the above
    make check

Running unit tests and generating a coverage report.

    pytest tests/ 
    coverage-badge -f -o coverage.svg

Updating Packages
-----------------
    conda update --all
    conda env export > environment.yaml

Documentation
-------------
We use sphinx to document our pipeline and code.

    sphinx-quickstart

Docker Image Development
------------------------

To run the SEEK pipeline in Docker, first follow instructions to install `Docker <https://docs.docker.com/get-docker/>`_.

**NOTE: You will need approximately at least 8-9 GB free disc space to run the Docker container.**

To setup the container in your system, you will pull pre-built Docker images from
neuroseek's `Docker Hub <https://hub.docker.com/orgs/neuroseek/repositories>`_.

.. code-block:: bash

    $ make pull-all

This will now pull all Docker containers needed to run ``seek`` to your local machine.

Now if you type in ``docker container ls``\,
you should see the corresponding container.

.. code-block:: bash

   # turn recipe to image
   docker build <image_container_name>

   # turn image to containeer
   docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash

:doc:`To better understand how we use Docker, see our Docker playbook <docker_playbook>`
