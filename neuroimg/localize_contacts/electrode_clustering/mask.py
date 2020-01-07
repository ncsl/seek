import re

import nibabel as nb
from nibabel.affines import apply_affine
import numpy as np
import numpy.linalg as npl


class MaskVolume:
    @classmethod
    def apply_mask(self, brain_img, mask_img):
        """
        Applies a binarized brain mask to an input CT image file

        Parameters
        –---------
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
        brain_img_data = brain_img.get_fdata()
        mask_data = mask_img.get_fdata()

        mask_data[mask_data > 0] = 1
        masked_brain = np.multiply(brain_img_data, mask_data)

        masked_brain_img = nb.Nifti1Image(masked_brain, brain_img.affine)

        return masked_brain_img

    @classmethod
    def mask_electrodes(self, elec_coords_mm, brainmasked_ct_img):
        """
        Filters out electrodes that do not fall within brain
        matter of a CT image.

        Parameters
        –---------
            elec_coords_mm: dict(str: ndarray)
                Dictionary that maps contact label to corresponding
                xyz coordinates.

            brainmasked_ct_img: NiBabel image object
                A brain-masked CT image in the form of NiBabel image object.

        Returns
        -------
            elecvoxels_in_brain: dict(str: ndarray)
                A dictionary of contact coordinates in CT voxels that fall
                within the binarized brain mask. The keys are individual
                contact labels, the values are the corresponding coordinates
                in CT space.
        """
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        brainmasked_ct_affine = brainmasked_ct_img.affine
        inv_affine = npl.inv(
            brainmasked_ct_affine
        )  # Obtain inverse affine matrix to transform from xyz to CT voxel

        # Filter out channels (in CT space) not within brain mask
        elecvoxels_in_brain = {}
        for label, contact in elec_coords_mm.items():
            vox = apply_affine(inv_affine, contact)
            idx = tuple(map(int, vox))
            if brainmasked_ct_data[idx]:
                elecvoxels_in_brain[label] = vox

        return elecvoxels_in_brain

    @classmethod
    def group_contacts(self, elec_in_brain):
        """
        Groups individual contacts by the electrode to which they correspond.
        Sorts the contacts using the corresponding labels.

        Parameters
        –---------
            elec_in_brain: dict(str: ndarray)
                Dictionary of contact coordinates in CT voxels that fall within
                the brain matter.

        Returns
        -------
            labeled_contacts: dict(str: dict(str: ndarray))
                Dictionary of contacts grouped by electrode. An electrode name
                maps to a dictionary of contact labels and corresponding
                coordinates. The dictionary is in sorted order based on these
                labels.
        """
        labeled_contacts = {}

        for label, coord in elec_in_brain.items():
            elecname = re.findall(r"[A-Za-z']+", label)[0]
            labeled_contacts.setdefault(elecname, {})[label] = coord

        for elec in labeled_contacts:
            sorted_chans = sorted(
                labeled_contacts[elec].items(),
                key=lambda x: int(re.findall(r"\d+", x[0])[0])
            )
            labeled_contacts[elec] = dict(sorted_chans)

        return labeled_contacts
