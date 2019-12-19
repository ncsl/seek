# Localizing Electrodes SEEG and ECoG

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

