# Stereotactic ElectroEncephalography Kit (SEEK)

SEEK is designed for epilepsy researchers that would like to apply anatomical analyses to their
iEEG analysis of patients. Although it is built for epilepsy researchers, any researcher that has
a similar problem as described below can probably use it.

The SEEK workflow attempts to help you automate as much of the processing of raw 
neuroimaging data (pre-surgical T1 MRI, post-implant CT, and optionally post-surgical T1 MRI). 
The pre-surgical T1 MRI are used to reconstruct a structural parcellation of the subject's brain.
The post-implant CT are used to localize implanted electrode coordinates. The post-surgical T1 MRI 
is used to determine the regions that were surgically resected/ablated in a procedure.

## TODO

1. Host documentation pages for seek, visualization engine on netlifly or something. Have them link
to each other.
2. Is thalamic/hippocampus/brainstem segmentation before/after aseg2srf.sh call which generates ?
3. Run a full example through each.


## Data Organization
We attempt to adhere to the BIDS specification. Original raw source data is stored in 
``sourcedata/<subject_id>``. Transformations (``.xfm`` files), image files (``.nii``) are
stored in the BIDS root when applicable. FreeSurfer derivatives are stored in ``derivatives/freesurfer/<subject_id>``.

See https://bids-specification.readthedocs.io/en/stable/ for more information on 
the data organization.

## Data Input
All raw source data should be in either Dicom (``.dcm``), or Nifti (``.nii.gz``, ``.nii``) format 
to start out with. 

## Setup and Installation
Please refer to the repository ``README`` file for instructions on setting up 
your workstation environment to run the SEEK ``snakemake`` workflows.

# The Workflows
SEEK exposes a number of workflows for the user to run. We sequentially describe the 
recommended workflow for any user starting with the T1/CT raw image files (in ``.dcm``).

## Setting environment variable
First before proceeding you must set your environment variable ``SEEKHOME`` to be the
directory path that ``seek`` repository lives. For example:

    export SEEKHOME=/Users/adam2392/Documents/seek

In addition, SEEK will use `blender` in the visualization prep workflow. Please 
download blender and set the environment variable to its runtime path.



FreeSurfer uses multiple processors to run. One can set this using the environment variable.

    # Use max processing power on machine
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$(nproc)

## Setting configuration
Inside the root of the repository directory, you can set configuration variables to match your dataset. 
Please modify ``config/localconfig.yaml`` and ``config/subjects.tsv`` files. These work by defining a 
``bids_root`` where it will be the path to the BIDS dataset. We will assume that under
 ``bids_root/sourcedata/<subject_id>`` will be raw image files. 
 
 By setting the configuration files correctly, the workflows will assume the following data directory 
 structure:
 
    <bids_root_path>/
        sourcedata/
            <subject_id>/
                premri/
                postct/
                postmri/ (optional)
        derivatives/
            freesurfer/
                <subject_id>/
                    acpc/
                    elecs/
                    CT/
                    ...
        sub-<subject_id>/
            anat/
                ...
            ct/
                ...
            ieeg/
                ...
        ...
        
            
## 1. ``recon_workflow``
This workflow has input as the T1 MRI pre-surgical dicoms/nifti.

The first workflow is intended to i) prep the T1 MRI data and ii) run FreeSurfer commands to generate
the structural reconstruction process. 

Prepping the data is done in ``recon_workflow/rules/prep.smk``, where FSL's ``robustfov`` command
and ``acpcdetect`` are used to create a field of view (fov) that crops out the neck to allow ACPC 
(anterior commissure and posterior commissure) alignment using ``acpcdetect``. We use FSL 
to create a robust fov first because ``acpcdetect`` fails otherwise. Note that output files 
from this process are created in ``derivatives/freesurfer/<subject_id>/acpc`` and then BIDS-compliant
images will be saved in the corresponding BIDS directory.

The FreeSurfer reconstruction workflow will essentially wrap ``recon-all`` and some downstream
commands to get all your files as desired into the derivatives and BIDS root folders.

To run the workflow:

    cd workflow/recon_workflow/
    snakemake -n  # dry run
    snakemake --cores 1  # run the workflow

## 2. ``prep_localization``
This workflow's input is the post-implantation CT dicoms/nifti.

This will perform various `mrconvert` commands to get the image into the BIDS respective 
folders.

    cd prep_localization_workflow/
    snakemake -n  # dry run
    snakemake --cores 1  # run the workflow
    
## 2. Manual step: Contact localization on CT image using FieldTrip Toolbox
In this step, the preimplantation CT dicoms/nifti image is required and also some metadata. 
Users will be required to know the electrode labels of all the channels that are implanted, and 
**also be able to determine** left/right hemisphere on the CT image. One can do so by taking 
a look at clinical schematics of the implantation and knowing the distribution of electrodes on
either hemisphere of the brain. Thus, they can differentiate the left side from right side. 
This is important in the Fieldtrip GUI.

The users are recommended to run FieldTrip Toolbox's localization GUI, with a script inside
``workflow/scripts/electrode_localization.mat``. In order for this script to work, you will need 
to modify the ``bids_root`` path and also add to path the Fieldtrip Toolbox, which can be downloaded 
at https://www.fieldtriptoolbox.org/download/. We used ``fieldtrip-20190129`` in our tests and
development, but we don't foresee future versions impacting the workflow significantly.

Note that this script will require as input the paths to the BIDS T1 MRI, and CT image. 
It will then ask the user to:

1. Determine coordinate system of the anatomical T1 MRI: If you followed steps earlier to align 
the T1 MRI to ACPC, then it should be in RAS coordinates. You can double-check at this point.
2. Determine the coordinate system of the CT image and add landmarks: You will be 
asked to click on approximate nasion, LPA, RPA and z-point landmarks in the CT image to provide 
points of reference when performing coregistration. This step is **VERY IMPORTANT** because the LPA 
and RPA are generally only distinguishable by knowing the implantation schematic of the patient. 
Alternatively, one can hope that your CT data has a marker delineating left from right 
side of the brain. Note that this step is crucial because coregistration algorithms from ``SPM12``, 
or ``FSL flirt`` will work on coregistering the CT to the T1 MRI, but it cannot determine 
which side of the CT is left or right, so without accurately placing the LPA and RPA, the coordinates 
of your electrodes will be on the wrong hemisphere.

After these steps, the script will run a coregistration of the CT to T1 MRI using ``SPM12`` 
and then it will display both images overlaid in blue/red. If the alignment looks approximately 
good enough, then you can proceed to localizing the electrode coordinates in the coregistered CT space.
The output of this will be xyz millimeter coordinates of each electrode in the anatomical MRI space.

*Reference:* For more information on Fieldtrip Toolbox and the GUI, see: https://www.fieldtriptoolbox.org/tutorial/human_ecog/#anatomical-workflow

### 2b. Advanced Usage - Semi-automated Electrode Localization

TODO. IN TESTING PHASE.

SEEK will provide an algorithm that attempts to localize all contacts 
on sEEG electrodes. In order to do so, one must localize at least 2 contacts 
per sEEG electrode in the above step. Then save the output, and run 
the corresponding algorithm script.

One can then open up the matlab script again to check the work of all the 
electrode locations.

Reference:
- Huynh C., Li A., et al., “Towards Automatic Localization and Anatomical 
Labeling of Intracranial Depth Electrodes in Brain Images”. 
IEEE Engineering in Medicine and Biology Conference, Montreal, Canada (2020).

## 3. ``coregistration_and_viz_workflow``
The input of this workflow is the coregistered CT image and electrode coordinates (in ``*_electrodes.tsv`` file 
format). 

By running this workflow, now many images will be mapped to the FreeSurfer ``T1.mgz`` 
file, which is in FreeSurfer MRI space. In addition, many transformation files will 
be saved that can be used in your analyses. They are self-described by their BIDS naming
scheme in the corresponding folders (``<subject_id>/ct/<fname>``). 

@CHRISTOPHER TO FILL IN DETAILS OF WHAT THIS IS DOING...
In addition, this workflow will generate Blender objects of the 

# Visualization Engine
Now that you have successfully ran through the SEEK workflow for your data, you will have 
created surface meshes that are anatomically labeled of the patient's brain, and also have
corresponding xyz (mm) locations and anatomical labels for each electrode in the same MRI space.
You can run a web-server now to visualize your data.

#### INSERT SCREENSHOT OF END

#### INSERT LINK TO A VERY SHORT VIDEO PERHAPS THAT CAN BE STORED ON GITHUB? CHRISTOPHER HAS ONE 

To run the visualization engine using Docker:
### INSERT MAKE COMMAND FOR RUNNING THE FINAL VISUALIZATION ENGINE 
    
    make 

# FAQ:
1. What if I want to map electrode coordinates to MNI space (e.g. and compare multiple subjects 
on the ``fsaverage`` subject space)?
    
    Answer: You can use the ``freesurfer/<subject_id>/mri/transforms/tal.xfm`` file and apply the 
    transformation to map your electrode coordinates mm (in subject's FreeSurfer MRI space) to MNI space.

2. 