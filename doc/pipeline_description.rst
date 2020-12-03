.. _PipelineDescription:

=========================
SEEK PIPELINE DESCRIPTION
=========================

The Snakemake rule-based workflow essentially abstracts the following rules and workflows. They are described in more
detail here:

#. (Optional) ACPC Alignment (pipeline/01-prep):

    Anterior commissure (AC) and posterior commissure (PC) alignment is generally a manual process.
    The midpoint of the anterior commissure (AC) is located first, serving as the origin of Talairach space.
    The brain is then rotated around the new origin (AC) so that the posterior commissure (PC) appears in the
    same axial plane as the anterior commissure. The connection of AC and PC in the middle of the brain forms the y-axis of
    the Talairach coordinate system. The x-axis runs from the left to the right hemisphere through AC.
    The z-axis runs from the inferior part of the brain to the superior part through AC.
    In order to further specify the x and z-axes, a y-z plane is rotated around the y (AC-PC) axis
    until it separates the left and right hemisphere (mid-sagittal plane).
    This completes the first part of Talairach transformation, often called AC-PC transformation.
    The obtained AC-PC space is attractive for individual clinical applications, especially presurgical
    mapping and neuronaviagation since it keeps the original size of the brain intact while providing a
    common orientation for each brain.
    
    .. code-block::
    
       acpcdetect -v -center-AC -output-orient LIP -no-tilt-correction -i ./T1.nii 

#. Reconstruction and Segmentation (pipeline/02-reconstruction):

    Reconstruction is the process of taking the raw images and matching the image volume with a specified atlas and then
    segmenting the brain into specific atlas regions. This can be done with FreeSurfer's recon-all command. 
    Note you want to convert things into Nifi format first for FS. 
    
    .. code-block::
    
       recon-all -i <patid_mriimg>.nii.gz -subjid <patid> -all

#. Coregistration (pipeline/03-coregistration)

    Coregistration maps different image volumes (e.g. T1 premri to T1 postop MRI) onto a 
    reference image volume. This puts the images into the same coordinate space, and also maps 
    voxels to voxels. This is generally done via some affine transformation computed 
    with external software packages.
    
    This can be done with FSL flirt command, or NDReg.    
    
    .. code-block::
    
       flirt -in {input.CT_NIFTI_IMG} \
                           -ref {input.PREMRI_NIFTI_IMG} \
                           -omat {output.output_affine_file} \
                           -out {output.CT_IN_T1PRE_NIFTI_IMG}

#. Contact Localization (pipeline/04-contact_localization: 

    Determines the locations (xyz) of the electrode channels in T1 MRI space (pipeline/contact_localization).
    This requires the user to first have preprocessed the CT scans (and optionally the T1 MRI). 

   .. code-block::

       matlab ./pipeline/contact_localization/matlab/run_localization_fieldtrip.m


    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 

#. Visualization (pipeline/05-visualization):

    This pipeline will visualize on a local web server the data.



Running Individual Workflows on Docker
--------------------------------------

Setup the sourcedata directory as follows:

.. code-block::

   sourcedata/{subject}/
                       - premri/*.dcm
                       - posmri/*.dcm
                       - postct/*.dcm


Edit the ``seek/pipeline/config/localconfig.yml`` file to include a list of
subject ids that you want to analyze.

Then run the following commands (assuming you built the container w/ Docker already):


#. Run reconstruction container:

   ..

      docker-compose run reconstruction /bin/bash


#. Run a rule, like "Prepping" the data

   Prep the data

   ..

      snakemake --snakefile ./pipeline/01-prep/Snakefile --cores 2


Workflow Steps Explained
------------------------

#. Setup data directory

   Setup your raw data directory for your patient to be read in:
    
    :: 

         study_name/
            /raw/
                /<patient_id>/
                    /premri/
                    /postct/
                    /postmri/ (optional)
     
    For more details, see BIDS: https://bids.neuroimaging.io/

#. Patient in FreeSurfer format explained:

    Here is how an example output `derivatives` directory would look like.

    ::
        derivatives/
            freesurfer/
                patient_id/
                    /mri/
                    /surf/
                    /label/
                    /stats/
                    /elecs/
                    /CT/
                    /ascii/
                    /Meshes/
                    /acpc/

    <patient_id> = The subject directory for data ran through FS (e.g. "umf001")
    /mri/ = Includes the mri-derived image transformations, including the original mri image volume.
    /surf/ = Includes the computed surface files for each hemisphere (rh and lh), such as white matter (wm), volume, thickness, pial, and smoothed surfaces
    /label/ = Includes derived labels for each surface mesh.
    /stats/ = Includes statistics computed for example for white matter, cortical volume.
    Additional Dirs Made Within to be compatible with FS
    /elecs/ = Localized contacts with xyz coordinates, anatomical mapping, etc.
    /CT/ = a directory to store the CT image volume and any transformations (e.g. mapped into T1 image volume)
    /ascii/ = ascii type files that show the subcortical volume.
    /Meshes/ = .mat files for the hemispheres and the triangular/vertices files for cortical and subcortical.
    /acpc/ = Anterior-commissure & posterior-commissure aligned image volumes. This is generally a common preprocessing step in many pipelines.


#. Change config.yaml file

    For the respective data directories of your data. This is under pipeline/config/localconfig.yaml

    * define `bids_root` directory

#. Run dry-run snakemake to make sure DAG job is constructed properly.

    Note, that you can only run snakemake commands after installing SnakeMake.

     .. code-block::

          snakemake -n # dry run
          snakemake # real run

#. Reconstruction

   .. code-block::

       cd pipeline/02-reconstruction
       snakemake -n
       snakemake

#. Coregistration

   .. code-block::

       cd pipeline/03-coregistration
       snakemake -n
       snakemake

#. Contact Localization

Note first, one should follow :ref:`LocalizationGuide`_ before running this.

.. code-block::

       cd pipeline/04-contact_localization
       snakemake -n
       snakemake

Snakemake Rules
---------------
Each of these workflows are enabled by a set of ``snakemake`` rules.
For an in-depth explanation on each particular rule, see `rules document <rules>`_.

Docker Usage in SEEK
--------------------
To heavily utilize Freesurfer, FSL, MRTrix3, and more, we make use of Docker.

:doc: `To better understand how we use Docker, see our Docker playbook <docker_playbook>.`
