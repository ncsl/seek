
PIPELINE DESCRIPTION
--------------------

The Snakemake rule-based workflow essentially abstracts the following rules and workflows.


#. 
   (Optional) ACPC Alignment:

   .. code-block::

       acpcdetect -v -center-AC -output-orient LIP -no-tilt-correction -i ./T1.nii 

#. 
   Reconstruction and Segmentation (pipeline/reconstruction):

    MRConvert to NIFTI format:

   .. code-block::

       mrconvert <mri_dir> <mri>.nii.gz


    Reconstruction is the process of taking the raw images and matching the image volume with a specified atlas and then
    segmenting the brain into specific atlas regions. 

    This can be done with FreeSurfer's recon-all command, or MRICloud. Note you want to convert things into Nifi format first for FS. 

   .. code-block::

       recon-all -i <patid_mriimg>.nii.gz -subjid <patid> -all

#. 
   Coregistration (pipeline/coregistration)

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

#. 
   (Manual Step): Determines the locations (xyz) of the electrode channels in T1 MRI space (pipeline/contact_localization).

    This requires the user to first have preprocessed the CT scans (and optionally the T1 MRI). 

   .. code-block::

       matlab ./pipeline/contact_localization/matlab/run_localization_fieldtrip.m


    This will run an ~10-15 minute process to have users determine how to localize the channels. Note that
    you will need the corresponding implantation map (i.e. PPT, some image drawn up by clinician, or the implantation knowledge).
    Deep channels (i.e. A1, B1, B'1, etc.) are usually in the brain, while the last channels of
    an electrode are near the skull. 

#. 
   Contributing Additional Transformations / Pipelines

    For developers contribbuting new rules into the pipeline, add explicit names for data 
    you want merged into the final output, so that new users can easily step in and understand
    what that data point is used for.

#. 
   Postprocess Analysis (pipeline/postprocess_analysis):

    TBD

#. 
   Resection/Ablation Estimation Volume (pipeline/resection_ablation_estimation):

    TBD

Running Your Own Image Patients
-------------------------------


#. 
   Setup your raw data directory for your patient to be read in:


   * 
     study_name/


     * /raw/

       * /patient_id/ 

         * /premri/
         * /postct/
         * /acpc/ (optional)

     For more details, see BIDS: https://bids.neuroimaging.io/

#. 
   Patient in FreeSurfer format explained:


   * derivatives/

     * freesurfer/

       * patient_id/ = The subject directory for data ran through FS (e.g. "umf001")

         * /mri/ = Includes the mri-derived image transformations, including the original mri image volume.
         * /surf/ = Includes the computed surface files for each hemisphere (rh and lh), such as white matter (wm), volume, thickness, pial, and smoothed surfaces
         * /label/ = Includes derived labels for each surface mesh. 
         * /stats/ = Includes statistics computed for example for white matter, cortical volume.
         * Additional Dirs Made Within to be compatible with FS
         * /elecs/ = Localized contacts with xyz coordinates, anatomical mapping, etc.
         * /CT/ = a directory to store the CT image volume and any transformations (e.g. mapped into T1 image volume)
         * /ascii/ = ascii type files that show the subcortical volume.
         * /Meshes/ = .mat files for the hemispheres and the triangular/vertices files for cortical and subcortical.
         * /acpc/ = Anterior-commissure & posterior-commissure aligned image volumes. This is generally a common preprocessing step in many pipelines.
         * /connectome/ = Any sort of connectome related files. For example, structural connectivity matrices used for The Virtual Brain.

#. 
   Change config.yaml file (local or cluster) to the respective
   data directories of your data. This is under pipeline/config/localconfig.yaml


   * define rawdata dir 
   * define FS output data dir (i.e. the FS_SUBJDIR)

#. 
   Run dry-run snakemake to make sure DAG job is constructed properly. Note, that you can only run snakemake commands after installing SnakeMake.


   * 
     you can run this in each of the subdirectories of pipeline/

     .. code-block::

          snakemake -n # dry run
          snakemake # real run

#. 
   Reconstruction

   .. code-block::

       cd pipeline/reconstruction
       snakemake -n
       snakemake

#. 
   Coregistration

   .. code-block::

       cd pipeline/coregistration
       snakemake -n
       snakemake

#. 
   Contact Localization

.. code-block::

       cd pipeline/contact_localization/matlab
       matlab
       <open run_localization_fieldtrip_v3.m>
       <change directories and variables>
       <run GUI>        
