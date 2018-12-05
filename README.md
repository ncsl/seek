# Neuroimaging Pipeline
By: Adam Li
10/4/18

# Setup and Installation
There are a couple of tools that you need to install in your system before everything is working. Preferably this is done via Docker, or Singularity, but if not, then:

1. Reconstruction
    * Freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
    * MRICloud (just send images here) (https://mricloud.org/)
2. Coregistration
    * NDReg (https://github.com/neurodata/ndreg)
    * FSL Flirt (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/)
3. Utility
    * MRTrix3 (https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html)

4. SPM, FieldTripToolbox, MATLAB


5. Conda env


    conda create -n <envname>
    conda install --file environment.yaml
    conda activate <envname>
    
    snakemake -n    


## Docker and Singularity
1. Freesurfer with FSL
2. MRTrix3
3. NDReG

TODO: Make sure SPM, FieldTripToolbox are imported as well
TBD

# Pipeline Description
At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

1. Reconstruction

2. Coregistration

3. Merging

4. (Manual Step): Determines the locations (xyz) of the electrode channels in T1 MRI space.
    
    This requires the user to first have preprocessed the T1 MRI, and the CT scans. Then
    the user must run coregistration, of which there are many options. Remember that coregistration
    is applying some deformation transformation on your input image (CT), such that it closely matches
    your reference image (T1 MRI).
    
    - NDREG
    - Flirt
    
    Once, coregistration is complete, the user should run:
           
        $ matlab ./matlab/run_localization_fieldtrip.m
        
    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 



### Pipeline Process Visualized
[DAG of Pipeline in Snakemake](./dag_preformat_pipeline.pdf)
