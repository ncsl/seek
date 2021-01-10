:orphan:

.. _installation:

INSTALLATION GUIDE
==================
``seek`` uses open-source third-party software to run the various workflows (e.g. `Freesurfer`_).
``seek`` itself is a wrapper using snakemake_. The best way to install the 3rd party software for ``seek`` usage
is via a Docker installation.

To fully install SEEK and run workflows, one will need to:

#. install SEEK repository
#. pull Docker containers

We outline some of these steps below. After you have set up everything (don't forget to
format your data repository according to our necessary format), then you can easily run
the snakemake workflows. For more information on running workflows after
installation, see :doc:`usage instructions <use>`.

Prerequisites
-------------
Install the following if you do not already have it:

    * `Docker`_
    * `Singularity`_ v3+ (you will need to install `Go`_ as well here)

If you decide to install things manually (without containers), we unfortunately
cannot guarantee things will work.

seek Installation
-----------------
There are a few ways to install seek itself.

.. code-block:: bash

    # clone repository locally
    $ git clone https://github.com/ncsl/seek
    $ python3.8 -m venv .venv
    $ pipenv install

Now to test that your installation worked, you can perform a dry-run of the first workflow.
Say you only have T1 `.dicom` images so far collected and formatted into a structured directory
as such:

.. code-block::

   {bids_root}/
        /sourcedata/
            /{subject}/
                - premri/*.dcm
                - postmri/*.dcm
                - postct/*.dcm

Next, you would want to modify the local configuration file as such to specify
different subjects you would like to run.

.. code-block::

    # modify config yaml
    $ cd seek/
    $ vim config/localconfig.yaml

    # modify the subjects table
    $ vim config/subjects.tsv

Finally, you can do a dry-run of the `recon_workflow <https://github.com/ncsl/seek/tree/master/workflow/recon_workflow>`_
by using snakemake.

.. code-block::

    $ cd workflow/recon_workflow/
    $ snakemake -n

You can also create a set of DAGS for yourself to visualize locally by running from the ``seek`` root
directory.

.. code-block::

    $ snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./recon_workflow.pdf;

This will in turn generate a ``.pdf`` file, which has a visual DAG of all the rules that will be run for the
subjects you specified in the ``subjects.tsv`` file.

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

Manual Installation (Not Recommended; See Docker)
-------------------------------------------------

For purposes of documentation and transparency to users, we outline here the manual installation process SEEK can take.
To install the SEEK pipeline manually, one must install the necessary python runtimes, as well as the necessary 3rd party
software.

Python Installations
^^^^^^^^^^^^^^^^^^^^

There are a couple of tools that you need to install in your system before everything is working. You ar recommended to use a Linux based OS. 
Follow links and tutorials on each respective tool to install. Preferably this is done via Docker, or Singularity, but if not, then:

Anaconda and Python3.6+: Conda (https://docs.anaconda.com/anaconda/install/)

This is mainly necessary to run snakemake, and any Python wrapper code.

.. code-block::

    conda env create -f environment.yml --name=seek
    source activate seek
    conda install sphinx sphinx-gallery sphinx_bootstrap_theme numpydoc black pytest pytest-cov coverage codespell pydocstyle
    pip install coverage-badge anybadge
    # dev versions of mne-python, mne-bids
    pip install --upgrade --no-deps https://api.github.com/repos/mne-tools/mne-python/zipball/master
    pip install --upgrade https://api.github.com/repos/mne-tools/mne-bids/zipball/master


Pipeline Installations (3rd Party Modules to Install)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Octave

    Runs open-source MATLAB-like functions. This runs various scripts for converting output files to object files for rendering visualizations.
    Follow: https://www.gnu.org/software/octave/#install

    .. code-block::

       brew install octave

#. Gawk_

    Runs command line tools.

#. Blender_

    Allows nice 3D mesh creations

#. Reconstruction (Freesurfer_)

    This step is necessary to generate a parcellation and surface reconstruction of the patient's brain.
    The general requirements is just a Linux, or OSX computer with enough RAM.
    Currently, this repo is designed to work with FreeSurfer.

#. Coregistration (`FSL Flirt`_)

    This step is necessary to map different imaging sessions together. Specifically, for this pipeline, we need it to map CT images to T1 MRI
    Note that as of 2019, installation still requires Python2, which should come in any Linux distribution.

     .. code-block::

          python2 <run_installer>

#. Utility (MRTrix3_)

#. SPM_ (preferably 12):

#. Contact-Localization Software (FieldTripToolbox, Img_Pipe, MATLAB)

   * FieldTripToolbox_

#. `ACPC Auto Detection (V2) <https://www.nitrc.org/projects/art/>`:


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