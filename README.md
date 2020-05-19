SEEK Pipeline (Stereotactic ElectroEncephalography Kit)
----------------------------------------------

[![CircleCI](https://circleci.com/gh/ncsl/seek.svg?style=svg)](https://circleci.com/gh/ncsl/seek)
[![Build Status](https://travis-ci.com/ncsl/seek.svg?token=6sshyCajdyLy6EhT8YAq&branch=master)](https://travis-ci.com/ncsl/seek)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
![GitHub](https://img.shields.io/github/license/ncsl/seek)
![GitHub last commit](https://img.shields.io/github/last-commit/ncsl/seek)
<a href="https://codeclimate.com/github/ncsl/seek/maintainability"><img src="https://api.codeclimate.com/v1/badges/2c7d5910e89350b967c8/maintainability" /></a>
![GitHub repo size](https://img.shields.io/github/repo-size/ncsl/seek)
[![DOI](https://zenodo.org/badge/160566959.svg)](https://zenodo.org/badge/latestdoi/160566959)
[![](https://images.microbadger.com/badges/version/neuroseek/seek.svg)](https://microbadger.com/images/neuroseek/seek "Get your own version badge on microbadger.com")

This repo describes Sarma/Crone lab effort to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT,
and iEEG data (ECoG, or SEEG). 

For incorporation of DTI data, see ndmeg: https://github.com/neurodata/ndmg

<!-- MarkdownTOC -->

- Features
- Setup and Installation
    - Modules to Install
- Data Organization
- Running Your Own Image Patients
- Pipeline Description
- Localizing Electrodes Process
    - Running Localization GUI
- Documentation
        - Pipeline Process Visualized
- References:

<!-- /MarkdownTOC -->

Features
--------
- [ ] Add support for MRICloud running using R-script. Possibly convert to Python script.
- [ ] Create unit and integration tests using pytest that test: pipeline in both snakemake and Python

Setup and Installation
--------
See [INSTALLATION GUIDE](doc/INSTALLATION.md)

### DOCKER

Setup: Note that the docker container names are:

    - neuroimg

To setup the container in your system:

    docker-compose up --build
    
    # run the container
    docker-compose up 
    
Now if you type in `docker container ls`, you should see the corresponding container.
    
    # turn recipe to image
    docker build <image_container_name>
    
    # turn image to containeer
    docker run -v $PWD/Data:/data -it -e bids_root=/data -e derivatives_output_dir=/data/derivatives --rm neuroimg_pipeline_reconstruction bash

# Creating persistent volumes

#### Reconstruction
1. Ensure that the data folder is set up as follows:
    - Data/sourcedata/neuroimaging/{subject}/
        - premri/*.dcm
        - posmri/*.dcm
        - postct/*.dcm
2. Build images:
    >
    <code>host:~# docker-compose up --build</code>
3. Run reconstruction container:
    >
    <code>host:~# docker-compose run reconstruction /bin/bash</code>
4. Prep the data
    > 
    <code>container:/neuroimg# snakemake --snakefile ./pipeline/01-prep/Snakefile --cores 2</code>
5.  Perform reconstruction
    >
    <code>container:/neuroimg# snakemake --snakefile ./pipeline/02-reconstruction/Snakefile --cores 2</code>
6. Perform coregistration
    >
    <code>container:/neuroimg# snakemake --snakefile ./pipeline/03-coregistration/Snakefile --cores 2</code>

#### Electrode localization (Bioimage suite)

1. Run localization container:
    >
    <code>host:~# docker-compose run localization ./start_bioimagesuite</code>

Data Organization
--------

We use BIDS. 
See https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy

Pipeline Description
--------
At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, 
regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

See [PIPELINE GUIDE](doc/PIPELINE_DESCRIPTION.md)

Semi-Automated Localizing Electrodes Process
-------- 
Localizing SEEG electrodes requires at least two contacts on each electrode to initialize the algorithm.
These can be say the deepest 2 contacts, or the entry point and target point (e.g. first and last contact on the electrode).

For ECoG data, we do not explicitly have a process outlined, but these are significantly easier since grids can
be easily interpolated.

See [LOCALIZATION_GUIDE](doc/LOCALIZATION_GUIDE.md)

Documentation and Testing
--------

See [Testing Guide](doc/TESTING_SETUP.md)
    
### Pipeline Process Visualized
[DAG of Pipeline in Snakemake](seek/neuroimg/pipeline/dag_neuroimaging_pipeline_reconstruction.pdf)

References:
--------
1. Recon-all. FreeSurfer. https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all#References
2. FSL Flirt. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT
3. MRTrix3. http://www.mrtrix.org/
4. Img_pipe. https://github.com/ChangLabUcsf/img_pipe
5. MRICloud. https://mricloud.org/
6. Snakemake. https://snakemake.readthedocs.io/en/stable/
7. FieldTrip Toolbox. http://www.fieldtriptoolbox.org/tutorial/human_ecog/


