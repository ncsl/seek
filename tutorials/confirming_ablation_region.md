# Confirming the Anatomical Region of the Ablation

This workflow is used to confirm that the annotated post MRI mask was mapped to the correct space by the following steps:

1. Manually annotate the ablation region on the post surgical MRI with Seg3D.
2. Use Slicer to orient the segmented MRI to freesurfer (FS) space.
3. Overlay the mask with the Desikan-Killany atlas and compare to clinical results.
4. Overlay the mask with electrode location and compare to clinical results. 

## 1. Manual Annotation

We use Seg3D (https://github.com/SCIInstitute/Seg3D/releases) to manually segment the postsurgical MRI, which creates a binary mask of the ablation region. 

Note that the postsurgical MRI should be in FS space. If not, map it to FS space using FLIRT. 


```python
flirt -in {ses-postsurgery_space-ACPC} -ref {ses-presurgery_space-fs} -out {ses-postsurgery_space-fs}
```

## 2. Mapping to FS Space

Because Seg3D uses its own coordinate system, we need to map the mask to FS space. This can be done using the 'Resample Image (BRAINS)' function in Slicer (https://download.slicer.org/).

Using the original postMRI (in FS space) as the reference image, we can map the mask to FS space. It may be a good idea to use an MRI viewer to ensure that the postMRI, the mask, and the reference atlas are all correctly aligned before moving on to the next step.

## 3. Overlay Mask with Atlas

Freesurfer provides three different types of atlases: Desikan-Killiany, DKT, and Destrieux. We use Desikan-Killany. More information on atlases can be found here: https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation.

Freesurfer also provides a lookup table (LUT) that assigns a number to a certain brain region. The full table is documented here:https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT.

Comparing the annotated voxels of the mask to the Desikan-Killany atlas, then finding the corresponding brain region in the Freesurfer LUT should give us a list of the brain regions that were ablated. This should match what was reported clinically.

## 4. Overlay Mask with Electrode Location

There are two ways we can overlay the mask with the electrodes:
1. Convert electrode locations to voxels
2. Convert mask to xyz coordinates

### Using Voxels

We can convert the electrode locations to voxels using the affine matrix of the presurgical MRI. Applying the inverse affine transformation then flooring the values will produce the voxels of where each electrode is located. 

Overlaying these numbers with the voxels of the mask produces a list of electrodes that were removed and the brain regions each electrode was associated with. This should match what was reported clinically. 

### Using XYZ Coordinates

We can convert the mask to xyz coordinates using the affine matrix of the postsurgical MRI. Applying the affine transformation will produce the xyz coordinates of the ablation region.

Overlaying these numbers with the xyz coordinates of the electrode locations produces a list of electrodes that were removed and the brain regions each electrode was associated with. This should match what was reported clinically.


```python

```
