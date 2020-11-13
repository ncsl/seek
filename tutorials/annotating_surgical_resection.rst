Guide to Segmenting and Annotating the Surgical Resected/Ablated Zone
---------------------------------------------------------------------

With `seek`, one can quickly get BIDs-compliant datasets that can then be used in downstream analysis.

Although we automate the coregistration mapping commands and file outputs between pre T1 and post T1 MRI
images, one still needs to create a mask of the post-surgical T1 MRI in order to annotate which voxels
were resected/ablated. Then by mapping this mask back into the pre-surgical T1 MRI space, then one can map
this into `FreeSurfer` space and get the voxels, brain regions and electrodes within the surgically
treated region.

