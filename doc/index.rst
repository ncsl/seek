What is SEEK?
=============

SEEK is a Snakemake-Python package that allows you to output
`BIDS <https://bids.neuroimaging.io/>`_\ -compatible datasets with the help of
Freesurfer_, FSL_, MRTrix3_, Blender_ and FieldTripToolbox_.

Why?
----
Using the above listed software to generate i) BIDS-compatible files and ii)
compute a pipeline that will obtain anatomical information for iEEG electrodes
is highly complex. There is a big learning curve, and custom bash scripts to generate
these files for many subjects are prone to error.

snakemake_ allows us to wrap all these rules using input and output files as a
checking system to create a Directed Acyclic Graph (DAG) of all the workflows connected.
All that is needed to start are sets of ``.dicom`` files for T1 MRI and CT (post iEEG
implantation). It also facilitates parallelization running of many subjects all at once.
If any part of the DAG fails for a specific subject, snakemake will be able to pick up where
it left off on the DAG.

What is the result?
-------------------

After running through the seek workflows, one will obtain BIDS-iEEG compatible ``*electrodes.tsv`` files
which contain anatomical information, and also numerous anatomical output files without
you worrying about learning the ins and outs of FreeSurfer, FSL and other software.

Finally, we pipe specific output to a visualization engine that is
housed in another `github repository <https://github.com/cronelab/ReconstructionVisualizer>`_.

Citing SEEK
-----------

If you use SEEK in your work, please cite our
`zenodo <https://zenodo.org/badge/latestdoi/160566959>`_.

Please also cite the following papers to credit BIDS and iEEG-BIDS:

- `BIDS <https://doi.org/10.1038/sdata.2016.44>`_
- `BIDS-iEEG <https://doi.org/10.1038/s41597-019-0105-7>`_


.. contents:: :local:
    :depth: 3

.. _Blender: https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/
.. _Freesurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
.. _FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/
.. _MRTrix3: https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html
.. _SPM: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
.. _FieldTripToolbox: http://www.fieldtriptoolbox.org/download/
.. _snakemake: https://snakemake.readthedocs.io/en/stable/