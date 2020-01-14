# Data Description
Here, I define the layout of data that will be optimal for running a patient. The minimal
amount of data needed is a T1 MRI scan (i.e. list of .dicom files).

CT and iEEG data files will be needed to perform electrode localization.

DTI .dicom files will be needed to perform tractography analysis.


#### Layout:
0. patient (raw)
   - premri
       - (list) dicoms
   - postct
       - (list) dicoms
   - dti
       - (list) dicoms

0b. freesurfer_output
   - patient (reconstruction)
       - (tree) files
   - outputfiles (stores all files we want to use in the next step)
   - success_file_flag.txt

1. patient (processed)
   - mri
       - nifti volume image: gives information about the T1 raw volume
   - ct
       - nifti volume image: gives information about the CT raw volume
   - dti (optional)
       - nifti volume image
   - elec
       - ct_in_t1_volume image: the CT volume coregistered to the T1 space
       - chans_voxels_xyz.txt: this is the list of voxel coordinates for each channel
       - chans_xyz.txt: this is the center of the estimated channel
   - surface (how to draw out the brain; note this is at the voxel level)
       - triangles.txt
       - vertices.txt
       - normals.txt
       - (TO ADD) voxels.txt 
   - parcellation (depends on atlas)
       - averaged:
           - region_mapping_cort_<atlas>.txt
           - region_mapping_subcort_<atlas>.txt
           - region_centres_<atlas>.txt
           - region_areas_<atlas>.txt
           - region_cortical_<atlas>.txt
           - region_average_orientations_<atlas>.txt
           - label_in_T1_<atlas> volume image (maps all voxels to an atlas region)
       - voxel (TO ADD):
           - region_mapping_cort_<atlas>.txt
           - region_mapping_subcort_<atlas>.txt
           - region_voxels_<atlas>.txt: the list of each voxel's label according to parcellation
           
   - connectome (depends on atlas)
       - tract_lengths_<atlas>.txt
       - weights_<atlas>.txt
       - gain_mat_<atlas>.txt (how to project regional activity -> sensor space | mainly for TVB)
   - seeg
       - edf
       - fif
   - scalp
       - edf
       - fif
   - clinical
       - ez_hypothesis_chans.txt
       - ez_hypothesis_<atlas>.txt
       - stores clinical metadata