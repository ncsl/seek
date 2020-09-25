from typing import List, Union

import numpy as np
import nibabel as nb


from mne_bids import BIDSPath
from mne_bids.tsv_handler import _from_tsv
from nibabel.affines import apply_affine
from nibabel.orientations import aff2axcodes


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
    """Map mm xyz coordinates into voxel space."""
    affine = img.affine
    inv_affine = np.linalg.inv(affine)
    img_data = img.get_fdata()

    # store surgical channels as a list
    surgical_chs = []

    # map coordinates to voxels and determine if they lay within segmentation volume
    for name, coord in zip(ch_names, ch_coords):
        new_coord = list(map(int, apply_affine(inv_affine, coord)))
        if img_data[new_coord] > 0:
            surgical_chs.append(name)
    return surgical_chs


def _compare_surgical_contacts(surg_chs, clinically_removed_contacts):
    not_found_chs = set(clinically_removed_contacts).difference(set(surg_chs))
    print("Difference: ", not_found_chs)
    return not_found_chs


if __name__ == "__main__":
    from pathlib import Path
    import nrrd

    bids_root = Path("/Users/adam2392/Dropbox/epilepsy_bids/")
    deriv_dir = Path(bids_root / "derivatives" / "freesurfer")
    subject = "la02"
    postsurg_dir = Path(deriv_dir / subject / "postsurgerymri")

    # get the segmented fpath
    segmented_fpath = Path(postsurg_dir / f"{subject}-surgical-segmentation.nrrd")

    # read the actual mask
    surgmask_arr, surgmask_hdr = nrrd.read(segmented_fpath)

    surgmask_arr = np.swapaxes(surgmask_arr, 1, 2)

    # get the original T1 img
    T1w_fpath = Path(postsurg_dir / "preT1.nii").as_posix()
    t1w_img = nb.load(T1w_fpath)
    t1w_affine = t1w_img.affine
    mm_to_voxel = np.linalg.inv(t1w_affine)
    ct_fpath = Path(deriv_dir / "CT" / "CT.nii")

    # get the post -> pre T1 affine
    post_to_pre_affine = np.loadtxt(Path(postsurg_dir / "fsl_postt1-to-t1_omat.txt"))
    pre_to_post_affine = np.linalg.inv(post_to_pre_affine)

    # print(surgmask_arr.shape, t1w_img.shape)
    # print(surgmask_hdr)
    # print(aff2axcodes(t1w_img.affine))
    # print(aff2axcodes(post_to_pre_affine))
    # print(nb.io_orientation(t1w_img.affine))

    # read in electrodes
    # subj_dir = make_bids_folders(
    #     bids_root=bids_root.as_posix(), subject=subject, session="veeg", make_dir=False
    # )
    # # get electrodes
    # electrodes_fpath = BIDSPath(
    #     subject=subject,
    #     session="veeg",
    #     processing="manual",
    #     acquisition="seeg",
    #     prefix=subj_dir,
    #     suffix="electrodes.tsv",
    # )
    # electrodes_tsv = _from_tsv(electrodes_fpath)

    ch_names = electrodes_tsv["name"]
    ch_coords = np.vstack(
        (electrodes_tsv["x"], electrodes_tsv["y"], electrodes_tsv["z"])
    ).T.astype(float)

    # convert CT coordinates to voxels

    # convert coordinates to voxels
    ch_voxs = np.array(
        [apply_affine(mm_to_voxel, coord) for coord in ch_coords]
    ).astype(int)
    # ch_voxs = np.swapaxes(ch_voxs)
    # ch_voxs = np.array([apply_affine(pre_to_post_affine, coord) for coord in ch_voxs]).astype(int)
    # print(ch_voxs)

    surgical_chs = []
    for i, (ch_name, ch_vox) in enumerate(zip(ch_names, ch_voxs)):
        # print(surgmask_arr[ch_vox[1], ch_vox[0], ch_vox[2]])
        if surgmask_arr[ch_vox[0], ch_vox[1], ch_vox[2]] == 1:
            surgical_chs.append(ch_name)

    print(ch_voxs)
    print(np.where(surgmask_arr > 0))

    print(surgical_chs)
    # print(ch_voxs)
