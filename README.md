# Neuroimaging Pipeline
By: Adam Li

Date: 10/4/18

This repo describes Sarma lab effort to pipeline explicitly a neuroimaging data workflow that involves T1 MRI, CT, DTI 
and iEEG data (ECoG, or SEEG). 

<!-- MarkdownTOC -->

- Setup and Installation
- Running Your Own Image Patients
    - Docker and Singularity
- Pipeline Description
- Localizing Electrodes Process
        - Pipeline Process Visualized

<!-- /MarkdownTOC -->

# Setup and Installation
There are a couple of tools that you need to install in your system before everything is working. Preferably this is done via Docker, or Singularity, but if not, then:

0. Anaconda and Python3.6+
    * Conda (https://conda.io/docs/user-guide/install/index.html)

1. Reconstruction
    * Freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
    * MRICloud (just send images here) (https://mricloud.org/)
    
2. Coregistration
    * NDReg (https://github.com/neurodata/ndreg)
    * FSL Flirt (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/)
    
3. Utility
    * MRTrix3 (https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html)

4. SPM, FieldTripToolbox, Img_Pipe MATLAB
    * SPM install (preferably 12): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
    * FieldTripToolbox (http://www.fieldtriptoolbox.org/download/)
    * Img_Pipe from the Chang Lab at UCSF will come as a submodule in this git repo. 
    
5. Conda env


    conda create -n <envname>
    conda activate <envname>
    conda env create -f environment_py3.yaml
    # conda install --file environment_py3.yaml
    cd pipeline/
    snakemake -n    

    # optionally separate install
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install numpy scipy pandas nibabel snakemake
    conda env export > environment_py3.yaml
    # conda env create -f environment_py3.yaml 

# Running Your Own Image Patients

1. Setup data directory for your patient in FreeSurfer format

2. Change config.yaml file (local or cluster) to the respective
data directories of your data.

3. Run "snakemake -n" to make sure DAG job is constructed propertly.

4. Run "snakemake"

## Docker and Singularity
1. Freesurfer with FSL
2. MRTrix3
3. NDReG

TODO: Make sure SPM, FieldTripToolbox are imported as well
TBD

# Pipeline Description
At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

1. Reconstruction and Segmentation:

    Reconstruction is the process of taking the raw images and matching the image volume with a specified atlas and then
    segmenting the brain into specific atlas regions. 
    
    This can be done with FreeSurfer's recon-all command, or MRICloud. 


2. Coregistration

    Coregistration maps different image volumes (e.g. T1 premri to T1 postop MRI) onto a 
    reference image volume. This puts the images into the same coordinate space, and also maps 
    voxels to voxels. This is generally done via some affine transformation computed 
    with external software packages.
    
    This can be done with FSL flirt command, or NDReg.    

3. Merging

    Freesurfer output and all other data is relatively complex. In order to make this repo usable
    by anyone with neuroscience knowledge, but lacking Reconstruction/Coregistration knowledge, a separate
    merging step will take place that extracts all relevant data and/or combines
    datapoints to allow for easy explicit I/O in matlab and python.
    
    For developers contribbuting new rules into the pipeline, add explicit names for data 
    you want merged into the final output, so that new users can easily step in and understand
    what that data point is used for.

4. (Manual Step): Determines the locations (xyz) of the electrode channels in T1 MRI space.
    
    This requires the user to first have preprocessed the T1 MRI, and the CT scans. Then
    the user must run coregistration, of which there are many options. Remember that coregistration
    is applying some deformation transformation on your input image (CT), such that it closely matches
    your reference image (T1 MRI).
    
    - NDREG
    - FSL
    
    Once, coregistration is complete, the user should run:
           
        $ matlab ./matlab/run_localization_fieldtrip_v3.m
        
    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 

# Localizing Electrodes Process

1. Use matlab script to get Voxel/MM coords in CT space

2. Apply coregistration transform matrix to coords to map to your MRI space.

3. (optional) Apply additional affine transformations to make sure your now T1 coordinates
are in an appropriate coordinate system and/or space, so that when using coords with atlas labels,
surface files and other T1-extracted image volumes, the coords are in the same 
language.

### Pipeline Process Visualized
[DAG of Pipeline in Snakemake](./pipeline/dag_neuroimaging_pipeline.pdf)
