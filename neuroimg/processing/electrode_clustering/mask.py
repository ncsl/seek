import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine


class MaskVolume():
    @classmethod
    def apply_mask(self, mask_data, brain_img_data):
        """
        Applies a binarized brain mask to an input CT image file
        :param mask_data: 3D array of brain mask
        :param brain_img_data: 3D array of CT brain scan
        :return: masked_brain: The masked CT image as a 3D array
        """

        mask_data[mask_data > 0] = 1
        masked_brain = np.multiply(brain_img_data, mask_data)

        return masked_brain

    @classmethod
    def filter_electrodes_bm(self, elec_coords_mm, brainmasked_ct_img):
        """
        Filters out electrodes that do not fall within brain matter of a CT image
        :param elec_coords_mm: A dictionary of contact coordinates in mm space
        :param brainmasked_ct_img: A brainmasked CT NiBabel image object
        :return: elecvoxels_in_brain: A dictionary of contact coordinates in CT voxels that fall
            within the binarized brain mask. The keys are individual contact labels, the values
            are the corresponding coordinates in CT space.
        """
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        brainmasked_ct_affine = brainmasked_ct_img.affine
        inv_affine = npl.inv(brainmasked_ct_affine)  # get the inverse affine matrix to go from xyz -> voxels

        # Convert contact xyz coordinates to CT voxels
        elec_coords_CTvox = {}
        for label, contact in elec_coords_mm.items():
            elec_coords_CTvox[label] = apply_affine(inv_affine, contact)

        # Filter out electrodes not within brain mask at the voxel level
        elecvoxels_in_brain = {}
        for label, contact in elec_coords_CTvox.items():
            if brainmasked_ct_data[int(contact[0]), int(contact[1]), int(contact[2])] != 0:
                elecvoxels_in_brain[label] = contact
        return elecvoxels_in_brain

    @classmethod
    def sort_contacts(self, elec_in_brain):
        """
        Groups the individual contacts by the electrode to which they correspond
        :param elec_in_brain: A dictionary of contact coordinates in CT voxels that fall
            within the brain matter.
        :return: voxels_per_electrode: A dictionary of contact coordinates in CT voxels that
            fall within the brain. The keys are electrode labels, and the values are lists of
            the coordinates of all the contacts that correspond to a given electrode.
        """
        voxels_per_electrode = {}
        for label, contact in elec_in_brain.items():
            if label[1] == "'":
                electrode_name = label[:2]  # If electrode name has a ' -> left hemisphere
            else:
                electrode_name = label[0]  # If electrode name has no '
            if electrode_name in voxels_per_electrode.keys():
                voxels_per_electrode[electrode_name].append(contact)
            else:
                voxels_per_electrode[electrode_name] = [contact]
