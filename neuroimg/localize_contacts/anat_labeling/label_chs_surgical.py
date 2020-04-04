from typing import List, Union

import numpy as np
import nibabel as nb

from nibabel.affines import apply_affine


def _apply_segmentation_mask(img, seg_mask_arr: np.ndarray) -> nb.Nifti2Image:
    img_data = img.get_fdata()

    # mask image array
    segmented_img_data = np.multiply(img_data, seg_mask_arr)

    # create Nifti image
    segmented_img = nb.Nifti2Image(segmented_img_data, img.affine, header=img.header)
    return segmented_img


def _get_surgical_contacts(
    img: nb.Nifti2Image, ch_names: List, ch_coords: Union[List, np.ndarray]
) -> List:
    affine = img.affine
    inv_affine = np.linalg.inv(affine)
    img_data = img.get_fdata()

    # store surgical channels as a list
    surgical_chs = []

    # map coordinates to voxels and determine if they lay within segmentation volume
    for name, coord in zip(ch_names, ch_coords):
        new_coord = list(map(int, apply_affine(affine, coord)))
        if img_data[new_coord] > 0:
            surgical_chs.append(name)
    return surgical_chs


def _compare_surgical_contacts(surg_chs, clinically_removed_contacts):
    not_found_chs = set(clinically_removed_contacts).difference(set(surg_chs))
    print("Difference: ", not_found_chs)
    return not_found_chs
