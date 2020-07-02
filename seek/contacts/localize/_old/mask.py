import nibabel as nb
import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine


class MaskVolume:
    """Class of functions of masking a neuroimaging volume."""

    @classmethod
    def apply_mask(self, brain_img: nb.Nifti2Image, mask_img: nb.Nifti2Image):
        """
        Apply a binarized brain mask to an input CT image file.

        Parameters
        ----------
            brain_img: NiBabel Image Object
                CT brain image - NIFTI formatting standard.

            mask_img: NiBabel Image Object
                Corresponding brain mask - NIFTI formatting standard.

        Returns
        -------
            masked_brain_img: NiBabel Image Object
                masked CT NiBabel Image Object in the same voxel space
                as input.
        """
        if brain_img.shape != mask_img.shape:  # pragma: no cover
            raise ValueError(
                "Brain image and Mask image shapes need to be the "
                f"same. You passed in {brain_img.shape} and  {mask_img.shape}."
            )

        # obtain image arrays
        brain_img_arr = brain_img.get_fdata()
        mask_arr = mask_img.get_fdata()

        # binarize mask to 0 and 1
        mask_arr[mask_arr > 0] = 1

        # multiply element wise brain array and mask array
        masked_brain = np.multiply(brain_img_arr, mask_arr)

        # return nibabel image
        masked_brain_img = nb.Nifti2Image(masked_brain, brain_img.affine)
        return masked_brain_img

    @classmethod
    def mask_electrodes(self, elec_coords_mm, brainmasked_ct_img):
        """
        Filter out electrodes that do not fall within brain matter of a CT image.

        Parameters
        ----------
            elec_coords_mm: dict(str: ndarray)
                Dictionary that maps contact label to corresponding
                xyz coordinates.

            brainmasked_ct_img: NiBabel image object
                A brain-masked CT image in the form of NiBabel image object.

        Returns
        -------
            elec_coords_mm: dict(str: ndarray)
                A dictionary of contact coordinates in CT voxels that fall
                within the binarized brain mask. The keys are individual
                contact labels, the values are the corresponding coordinates
                in CT space.
        """
        # brain masked data has 0 outside the brain
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()

        # Obtain inverse affine matrix to transform from xyz to CT voxel
        brainmasked_ct_affine = brainmasked_ct_img.affine
        inv_affine = npl.inv(brainmasked_ct_affine)

        # Filter out channels (in CT space) not within brain mask
        _elec_coords_mm = {}
        for label, ch_coord in elec_coords_mm.items():
            vox = apply_affine(inv_affine, ch_coord)
            idx = tuple(map(int, vox))
            if brainmasked_ct_data[idx] != 0:
                _elec_coords_mm[label] = ch_coord
        elec_coords_mm = _elec_coords_mm
        return elec_coords_mm
