# Localizing Electrodes SEEG

## User Experience Requirements
1. Voxel level analysis is required because
    - how many voxels correspond to a contact?
    - plotting voxel point clouds
2. Need voxel readable file
    - voxel txt file [x,y,z,<brain_id>]
    - brain_id Look up table file
    - store as a .mat and .npz file
3. Store voxel level files per atlas
4. Store 

	
3. store voxel level xyz coords of channels per atlas
- voxel_atlas.txt (<x,y,z>, label_index)
- label_index_atlas_lut.txt (label_index, atlas_brain_region)
- channel_voxel.txt (channels: <x,y,z>) [note; not atlas dependent]
- do surface files for voxel level -> look at aparc+aseg.mgz file
- Also, create the averaged files:
	- channel_xyz_avgpoint.txt (

## Implementation
1. Read in the aprac+aseg.mgz

## Tests
1. test failed data
	- Macauley/Kristin to insert fake data
2. test patient
	- LA02
	- EFRI03
3. test:
	- test plotting brain regions
	- test plotting seeg channels
	- test getting brain_region -> seeg channels


1. Input
- take CT Images
- take MRI Image, brain envelope 
- 

Methods:
- map brain envelope to the CT space and mask the skull out
- applying simple threshold 
	- can get clusters of bright spots that appear for varying different thresholds
	- 
- computing distance between each dot with all other clusters
	- distance matrix
	- use another threshold to decide which are connected and not
	- loop through and check if they are in a line or not
	- now each grouping is a vector of channels that represent an electrode
(note that it does not handle the case when the channels are not put into a line)


2. Output
- output list of channels grouped by electrodes
- excel sheet that can group these electrodes
	- missing electrode when algorithm fails
	-> can loop back to read into the function again with an entry point
	
## Data Description
Layout:
0. patient (raw)
   - premri
       - (list) dicoms
   - postct
       - (list) dicoms
   - dti
       - (list) dicoms

0b. freesurfer_output / reconstruction_output
   - patient (reconstruction)
       - (tree) files
   - tempdir (stores all files we want to use in the next step)
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