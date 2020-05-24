# INSTALLATION GUIDE
To install the SEEK pipeline, one must install the necessary python runtimes, as well as the necessary 3rd party 
softwares. 

The best way to do so is via a Docker installation.

<!-- MarkdownTOC -->

- Docker Installation
- Manual Installation
    - Python Installations
    - Pipeline Installations \(3rd Party Modules to Install\)

<!-- /MarkdownTOC -->

## Docker Installation
To run the SEEK pipeline in Docker, first follow instructions to install [Docker](https://docs.docker.com/get-docker/).

**NOTE: You will need approximately at least 8-9 GB free disc space to run the Docker container.**

To setup the container in your system:

    # build the composition in `docker-compose.yml`
    docker-compose up --build
    
    # run the container
    docker-compose up 
    
Now if you type in `docker container ls`, you should see the corresponding container.
    
    # turn recipe to image
    docker build <image_container_name>
    
    # turn image to containeer
    docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash

### Running Individual Workflows on Docker

#### Reconstruction Workflow
Setup the sourcedata directory as follows:

 **sourcedata/**
       
    /{subject}/
        - premri/*.dcm
        - posmri/*.dcm
        - postct/*.dcm
        
 
Edit the `seek/pipeline/config/localconfig.yml` file to include a list of 
subject ids that you want to analyze.

Then run the following commands (assuming you built the container w/ Docker already):

   1. Run reconstruction container:
   
        > docker-compose run reconstruction /bin/bash

   2. Prep the data
   
        >snakemake --snakefile ./pipeline/01-prep/Snakefile --cores 2
                              
   3.  Perform reconstruction
        
        > snakemake --snakefile ./pipeline/02-reconstruction/Snakefile --cores 2
   
   4. Perform coregistration
   
        > snakemake --snakefile ./pipeline/03-coregistration/Snakefile --cores 2

#### Electrode Localization Workflow
TBD

#### Visualization of Localized Electrodes
TBD



## Manual Installation
### Python Installations
There are a couple of tools that you need to install in your system before everything is working. You ar recommended to use a Linux based OS. 
Follow links and tutorials on each respective tool to install. Preferably this is done via Docker, or Singularity, but if not, then:

Anaconda and Python3.6+ :
   * Conda (https://docs.anaconda.com/anaconda/install/)
   * This is mainly necessary to run img_pipe (ECoG localization with Chang-Lab repo), snakemake, and any Python wrapper code
    
    i. 
        conda env create -f environment.yml --name=seek
        source activate seek
        conda install sphinx sphinx-gallery sphinx_bootstrap_theme numpydoc black pytest pytest-cov coverage codespell pydocstyle
        pip install coverage-badge anybadge
        # dev versions of mne-python, mne-bids
        pip install --upgrade --no-deps https://api.github.com/repos/mne-tools/mne-python/zipball/master
        pip install --upgrade https://api.github.com/repos/mne-tools/mne-bids/zipball/master
        
    ii. Conda env
    
        # create environment
        conda create -n seek
        conda activate seek
        
        # optionally separate install
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda install numpy scipy matplotlib scikit-learn scikit-image pandas seaborn nibabel mne snakemake mne-bids flask
        conda install pytest black check-manifest pytest-cov pydocstyle
        
        # check if installation worked
        cd neuroimg/pipeline/reconstruction/
        snakemake -n    
   
Pip and setup.py install
    
    # run installation via setup.py
    make inplace-all
    
    # install testing functionality
    make install-tests


### Pipeline Installations (3rd Party Modules to Install)

1. Octave

Runs open-source. This runs various scripts for converting output files to object files for rendering visualizations.
Follow: https://www.gnu.org/software/octave/#install
    
    brew install octave
       
2. Gawk
    Runs command line tools.
    https://brewinstall.org/Install-gawk-on-Mac-with-Brew/

3. Blender
    https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/

4. Reconstruction
    * Freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
    * This step is necessary to generate a parcellation and surface reconstruction of the patient's brain. The general requirements is just a 
    Linux, or OSX computer with enough RAM. Currently, this repo is designed to work with FreeSurfer.
    
5. Coregistration
    * FSL Flirt (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/)
    * This step is necessary to map different imaging sessions together. Specifically, for this pipeline, we need it to map CT images to T1 MRI
    * Note that as of 2019, installation still requires Python2, which should come in any Linux distribution.
        
            python2 <run_installer>
    
6. Utility
    * MRTrix3 (https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html)

7. SPM 
    * SPM install (preferably 12): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

8. Contact-Localization Software (FieldTripToolbox, Img_Pipe, MATLAB)
    * FieldTripToolbox (http://www.fieldtriptoolbox.org/download/)
    * Img_Pipe from the Chang-Lab at UCSF will come as a submodule in this git repo. This heavily handles ECoG data only.
 
9. ACPC Auto Detection (V2):
    * https://www.nitrc.org/projects/art/
    
10. (Optional) Cloud Reconstruction (MRICLOUD):
    * MRICloud (cloud based soln; just send images here) (https://mricloud.org/)
    * the nice thing is that this usually works even when FS fails (e.g. the T1 MRI image isn't good enough quality, or there is a major lesion, etc.).


11. (Optional) Nonlinear Registration NDREG:
    * NDReg (https://github.com/neurodata/ndreg)
