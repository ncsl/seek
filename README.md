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

# Features
- [x] Keep FS style directory, but instead add directories there if necessary. Remove fat code that simply copies and pastes data files
- [ ] Add Travis.CI testing
- [x] Convert pipeline into running on the FS directory structure
- [ ] Add support for MRICloud running using R-script. Possibly convert to Python script.
- [ ] Create Python pipeline for running everything, or submit PR to Chang Lab's to run SEEG.
- [ ] Create unit and integration tests using pytest that test: pipeline in both snakemake and Python
- [ ] Add visualization features for SEEG that can be synced into Chang's lab repo. 
- [ ] Add Docker container for project

# Setup and Installation
There are a couple of tools that you need to install in your system before everything is working. You ar recommended to use a Linux based OS. 
Follow links and tutorials on each respective tool to install. Preferably this is done via Docker, or Singularity, but if not, then:

0. Anaconda and Python3.6+ 
    * Conda (https://docs.anaconda.com/anaconda/install/)
    * This is mainly necessary to run img_pipe (ECoG localization with Chang Lab repo), snakemake, and any Python wrapper code
    
    i. Conda env
    
    
        # probably doesn't work
        conda create -n <envname>
        conda activate <envname>
        conda env create -f environment.yml
        # conda install --file environment.yml
        cd pipeline/
        snakemake -n    
    
        # optionally separate install
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda install numpy scipy pandas nibabel snakemake
        conda install -c conda-forge mne 
        conda install -c flyem-forge/label/upgrade201904 marching_cubes
        conda env export > environment.yml
        # conda env create -f environment_py3.yaml 

    
   ii. Update all submodules


        git submodule update --init --recursive 
        
        
1. Reconstruction
    * Freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
    * This step is necessary to generate a parcellation and surface reconstruction of the patient's brain. The general requirements is just a 
    Linux, or OSX computer with enough RAM. Currently, this repo is designed to work with FreeSurfer.
    
2. Coregistration
    * FSL Flirt (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/)
    * This step is necessary to map different imaging sessions together. Specifically, for this pipeline, we need it to map CT images to T1 MRI
    * Note that as of 2019, installation still requires Python2, which should come in any Linux distribution.
        
            python2 <run_installer>
    
3. Utility
    * MRTrix3 (https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html)

4. SPM 
    * SPM install (preferably 12): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

5. Contact-Localization Software (FieldTripToolbox, Img_Pipe, MATLAB)
    * FieldTripToolbox (http://www.fieldtriptoolbox.org/download/)
    * Img_Pipe from the Chang Lab at UCSF will come as a submodule in this git repo. This heavily handles ECoG data only.
 
6. (Optional) Cloud Reconstruction (MRICLOUD):
    * MRICloud (cloud based soln; just send images here) (https://mricloud.org/)
    * the nice thing is that this usually works even when FS fails (e.g. the T1 MRI image isn't good enough quality, or there is a major lesion, etc.).

7. (Optional) ACPC Auto Detection:
    * https://www.nitrc.org/projects/art/
    
8. (Optional) Nonlinear Registration NDREG:
    * NDReg (https://github.com/neurodata/ndreg)

# Running Your Own Image Patients

0. Setup your raw data directory for your patient to be read in:
    
    * /patient_id/ 
        - /premri/
        - /postct/
        - /dti/ (optional)

1. (Optional) Setup data directory for your patient in FreeSurfer format:

    * /patient_id = The subject directory for data ran through FS (e.g. "umf001")
        - /mri/ = Includes the mri-derived image transformations, including the original mri image volume.
        - /surf/ = Includes the computed surface files for each hemisphere (rh and lh), such as white matter (wm), volume, thickness, pial, and smoothed surfaces
        - /label/ = Includes derived labels for each surface mesh. 
        - /stats/ = Includes statistics computed for example for white matter, cortical volume.
        - Additional Dirs Made Within to be compatible with FS
        - /elecs/ = Localized contacts with xyz coordinates, anatomical mapping, etc.
        - /CT/ = a directory to store the CT image volume and any transformations (e.g. mapped into T1 image volume)
        - /ascii/ = ascii type files that show the subcortical volume.
        - /Meshes/ = .mat files for the hemispheres and the triangular/vertices files for cortical and subcortical.
        - /acpc/ = Anterior-commissure & posterior-commissure aligned image volumes. This is generally a common preprocessing step in many pipelines.
        - /connectome/ = Any sort of connectome related files. For example, structural connectivity matrices used for The Virtual Brain.
        
2. Change config.yaml file (local or cluster) to the respective
data directories of your data. This is under pipeline/config/localconfig.yaml
    
    * define rawdata dir 
    * define FS output data dir (i.e. the FS_SUBJDIR)

3. Run dry-run snakemake to make sure DAG job is constructed properly. Note, that you can only run snakemake commands after installing SnakeMake.
    
    * you can run this in each of the subdirectories of pipeline/
        
            snakemake -n # dry run
            snakemake # real run
            
4. Reconstruction

        cd pipeline/reconstruction
        snakemake -n
        snakemake

5. Coregistration

        cd pipeline/coregistration
        snakemake -n
        snakemake
        
6. Contact Localization

        
        cd pipeline/contact_localization/matlab
        matlab
        <open run_localization_fieldtrip_v3.m>
        <change directories and variables>
        <run GUI>        

## Docker and Singularity
1. Freesurfer with FSL
2. MRTrix3
3. NDReG

TODO: Make sure SPM, FieldTripToolbox are imported as well
TBD

# Pipeline Description
At a high level, this pipeline is taking neuroimaging data of a patient to produce usable data about the brain's geometry, regional parcellation into atlas regions, connectivity between brain regions measured by white matter tracts, and channel localization in MRI space.

0. (Optional) ACPC Alignment:
    
        acpcdetect -v -center-AC -output-orient LIP -no-tilt-correction -i ./T1.nii 

1. Reconstruction and Segmentation (pipeline/reconstruction):

    MRConvert to NIFTI format:
    
        mrconvert <mri_dir> <mri>.nii.gz
    
    Reconstruction is the process of taking the raw images and matching the image volume with a specified atlas and then
    segmenting the brain into specific atlas regions. 
    
    This can be done with FreeSurfer's recon-all command, or MRICloud. Note you want to convert things into Nifi format first for FS. 

        recon-all -i <patid_mriimg>.nii.gz -subjid <patid> -all
        
2. Coregistration (pipeline/coregistration)

    Coregistration maps different image volumes (e.g. T1 premri to T1 postop MRI) onto a 
    reference image volume. This puts the images into the same coordinate space, and also maps 
    voxels to voxels. This is generally done via some affine transformation computed 
    with external software packages.
    
    This can be done with FSL flirt command, or NDReg.    

        flirt -in {input.CT_NIFTI_IMG} \
                            -ref {input.PREMRI_NIFTI_IMG} \
                            -omat {output.output_affine_file} \
                            -out {output.CT_IN_T1PRE_NIFTI_IMG}

3. (Manual Step): Determines the locations (xyz) of the electrode channels in T1 MRI space (pipeline/contact_localization).
    
    This requires the user to first have preprocessed the T1 MRI, and the CT scans. Then
    the user must run coregistration, of which there are many options. Remember that coregistration
    is applying some deformation transformation on your input image (CT), such that it closely matches
    your reference image (T1 MRI). Once, coregistration is complete, the user should run:
           
        matlab ./matlab/run_localization_fieldtrip_v3.m
        
    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 

4. Contributing Additional Transformations / Pipelines

    For developers contribbuting new rules into the pipeline, add explicit names for data 
    you want merged into the final output, so that new users can easily step in and understand
    what that data point is used for.

5. Diffusion Tensor Imaging (pipeline/diffusion_tensor_imaging):

    TBD
    
6. Postprocess Analysis (pipeline/postprocess_analysis):

    TBD

7. Resection/Ablation Estimation Volume (pipeline/resection_ablation_estimation):

    TBD

# Localizing Electrodes Process

1. Use matlab script to get Voxel/MM coords in CT space

2. Apply coregistration transform matrix to coords to map to your MRI space.

3. (optional) Apply additional affine transformations to make sure your now T1 coordinates
are in an appropriate coordinate system and/or space, so that when using coords with atlas labels,
surface files and other T1-extracted image volumes, the coords are in the same 
language.

### Pipeline Process Visualized
[DAG of Pipeline in Snakemake](./pipeline/dag_neuroimaging_pipeline.pdf)

# References:
1. Recon-all. FreeSurfer. https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all#References
2. FSL Flirt. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT
3. MRTrix3. http://www.mrtrix.org/
4. Img_pipe. https://github.com/ChangLabUcsf/img_pipe
5. MRICloud. https://mricloud.org/
6. Snakemake. https://snakemake.readthedocs.io/en/stable/
7. FieldTrip Toolbox. http://www.fieldtriptoolbox.org/tutorial/human_ecog/