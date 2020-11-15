Guide to setup Testing Environment
==================================
Seek is a pipeline manager for the workflows involving SEEG electrode implantations,
T1 MRI and CT scans.

We use `snakemake` module to link the pipelines in a Directed Acyclic Graph (DAG).

Testing Data
------------

We utilize a testing dataset at https://github.com/adam2392/seek-testing-data to
perform unit and integration tests.

Testing Scipts
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