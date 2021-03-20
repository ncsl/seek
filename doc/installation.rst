:orphan:

.. _installation:

Installation Guide
==================
``seek`` uses open-source third-party software to run the various workflows (e.g. `Freesurfer`_).
``seek`` itself is a wrapper using snakemake_. The best way to install the 3rd party software for ``seek`` usage
is via a Docker installation.

To fully install SEEK and run workflows, one will need to:

#. install SEEK repository
#. pull Docker containers

We outline some of these steps below. After you have set up everything (don't forget to
format your data repository according to our necessary format), then you can easily run
the snakemake workflows. 

.. _data_organization:

Data Organization
-----------------

We assume data is organized in BIDS. If your data is not in BIDS, then ``seek-pipeline`` will not work for you.
See https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy

Before data is converted to BIDS in ``seek/pipeline/01-prep`` pipeline,
then ``sourcedata/`` should contain a semi-structured format of the neuroimaging data that will
be used in the ``snakemake`` workflows.

.. code-block::

    <bids root>/
        sourcedata/
            {subject}/
            - premri/*.dcm
            - posmri/*.dcm
            - postct/*.dcm

Prerequisites
-------------
Install the following if you do not already have it:

    * `Docker`_
    * `Singularity`_ v3+ (you will need to install `Go`_ as well here)

We recommend following the official documentation, but if you have trouble, 
we used a specific set of instructions that might work for you. See :ref:`singularity_help`.

If you decide to install things manually (without containers), we unfortunately
cannot guarantee things will work due to the fact that running pipelines on data 
require installing a variety of different softwares with different versions.

Installation: seek-pipeline
---------------------------

Once you have Docker_ and Singularity_ installed, you can install the ``seek-pipeline`` 
repository to run locally on your machine. There are a few ways to do so:

.. code-block:: bash

    # clone repository locally
    git clone https://github.com/ncsl/seek
    python3.8 -m venv .venv
    pipenv install

This will setup a Python3.8 virtual environment and install necessary dependencies
inside this virtual environment listed in the ``Pipfile``.

Setup configuration file
------------------------

Now to test that your installation worked, you can perform a dry-run of the first workflow.
Say you only have T1 `.dicom` images so far collected and formatted into a structured directory 
as specified in :ref:`data_organization`.

Next, you would want to modify the local configuration file as such to specify
different subjects you would like to run.

.. code-block::

    # modify config yaml
    cd seek/
    vim config/localconfig.yaml

    # modify the subjects table
    vim config/subjects.tsv

Dry-run and test configuration
------------------------------

Finally, you can do a dry-run of the `recon_workflow <https://github.com/ncsl/seek/tree/master/workflow/recon_workflow>`_
by using snakemake.

.. code-block::

    cd workflow/recon_workflow/
    snakemake -n

You can also create a set of DAGS for yourself to visualize locally by running from the ``seek`` root
directory.

.. code-block::

    snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./recon_workflow.pdf;

This will in turn generate a ``.pdf`` file, which has a visual DAG of all the rules that will be run for the
subjects you specified in the ``subjects.tsv`` file. An example would look like the following for one 
subject with identifier ``la02``:

.. image:: /_static/recon_workflow.png
    :width: 400
    :alt: Reconstruction workflow

Now that you have successfully installed ``seek-pipeline``, head on over to 
our :doc:`usage documentation <use>` to see how to run various workflows and 
add add additional configurations, such as parallelization.

.. _singularity_help:

Singularity Installation (for Linux)
------------------------------------
In order to run snakemake rules using Docker containers, you ``need`` Singularity.
Although installations differ and may evolve, we highlight an installation sequence
that worked for us locally. Otherwise, go to `Singularity`_'s installation page.

To install Singularity, we tested on version 3.7.0, but it should work
on any of the versions 3.5+.

When installing these, we used the Go version 1.15.6.
But minimally 1.13+ should work. Here are a few code snippets
for installing Go and then singularity.

.. code-block:: bash

    export VERSION=1.15.6 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

.. code-block:: bash

    echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

Now install singularity

.. code-block:: bash

    go get -d github.com/sylabs/singularity
    export VERSION=3.7.0 && # adjust this as necessary \
    mkdir -p $GOPATH/src/github.com/sylabs && \
    cd $GOPATH/src/github.com/sylabs && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd ./singularity && \
    ./mconfig

.. code-block:: bash

    ./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

.. _Gawk: https://brewinstall.org/Install-gawk-on-Mac-with-Brew/
.. _Blender: https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/
.. _Freesurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
.. _FSL Flirt: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/
.. _MRTrix3: https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html
.. _SPM: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
.. _FieldTripToolbox: http://www.fieldtriptoolbox.org/download/
.. _snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Docker: https://docs.docker.com/get-docker/
.. _Singularity: https://sylabs.io/guides/3.0/user-guide/installation.html
.. _Go: https://golang.org/doc/install