# Pipeline Readme

By: Adam Li

1/17/19


## Steps
1. Reconstruction
    
    Here, one performs reconstruction and segmentation using the T1 MRI scans of the full brain (i.e. preop).
    
    From here, one can either use FreeSurfer, or MRICloud. 
    
    Freesurfer will perform the entire pipeline using recon-all and create all files in their specified directory.
    The pipeline will then copy out files that are important for further downstream analysis. Things will include 
    segmentations based on Desikan-Killiany, and Destrieux atlas, 
    
2. Coregistration

    Here, one performs coregistration of different image volumes. Here the main image volumes that will
    be mapped are the CT -> T1 MRI scan. 
    
    Here, we implement rules to use Flirt FSL's algorithm, Robust Registration in FreeSurfer and Large Deformation 
    Diffeomorphic Metric Mapping Registration using ndreg.
    
    Under each rule, a mapped image volume and a txt file describing the affine transformation should be saved.
    
3. Contact Localization
    
    Here, one runs a semi-automated pipeline for localizing the contacts in CT space, running a clustering 
    algorithm and then mapping contacts into their corresponding MRI space.
    
    The pipeline will assume contacts are in native MRI space.

4. (Optional) Diffusion_Tensor_Imaging


5. (Optional) Resection_Ablation_Estimation

