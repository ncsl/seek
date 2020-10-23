Guide to Segmenting and Annotating the Surgical Resected/Ablated Zone
---------------------------------------------------------------------

With `seek`, one can quickly get BIDs-compliant datasets that can then be used in downstream analysis.

Although we automate the coregistration mapping commands and file outputs between pre T1 and post T1 MRI
images, one still needs to create a mask of the post-surgical T1 MRI in order to annotate which voxels
were resected/ablated. Then by mapping this mask back into the pre-surgical T1 MRI space, then one can map
this into `FreeSurfer` space and get the voxels, brain regions and electrodes within the surgically
treated region.

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