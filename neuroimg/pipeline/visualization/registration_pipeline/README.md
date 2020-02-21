# MRIs to Meshes

Converts a T1 MR into a segmented and parcellated 3D object file.

Each individual brain region is it's own mesh and can be queried/modified/disabled. Colormap is taking directly from Freesurfer's lookup table.

The gLTF binary file format allows for optimized viewing directly in browsers and VR/AR devices.
*Examples coming soon.*

This project is essentially feature complete.

### Dependencies:
- Freesurfer (90% of this project is just **recon-all**)
    - Which means it will take ~16ish hours to run.
- Various [Brainder](https://brainder.org/) scripts for mesh generation
    - Matlab/Octave
- Blender (Freesurfer OBJ to FBX/gLTF)
    - Comes with a python interpreter

- **Bioimage Suite**
    - This is where all of the coregistration/electrode localization is done. After warping the CT-marked electrodes to the Freesurfer-normalized T1.nii/.mgz export the BIS .mgrid as a .txt file and drop it in the *Freesurfer Subject Directory*/Subject/electrodes folder

### To do:
- Documentation/Instructions
- Investigate either embedding the former Unity/Holgram shader into the 3D scene on either the Blender side or ThreeJS side
- Option to just create a brain mesh and not do anything with electrodes
- Possibly remove the dependency on Blender
- Incorporate into the [neuroimg_pipe](https://github.com/adam2392/neuroimg_pipeline) repo
    - This will ideally replace the dependence on Bioimage suite (but it might be nice to keep that as an option)
- As long as the electrode file is property registered to the T1, any localizer should be able to be used. Could probably re-incorporate FieldTrip to make it more software-agnostic
- Clean up what actually gets exported. Some regions we probably don't care about.

### Instructions:

- Create a new folder in the Freesurfer subject directory and place the following files:
- T1.nii
- CT.nii
- electrodes.txt (from Bioimage Suite)
- T2.nii (optional)
- DTI and corresponding bvec/bval (optional)

### Issues:
- Freesurfer dev branch is required for hippocampal/amygdal/thalamic segmentation.
- Freesurfer 6.0 is required for trac-all
(No idea why. Something having to do with the bval/bvec file and/or flirt)

### Example:

![Example](/Picture.jpg)