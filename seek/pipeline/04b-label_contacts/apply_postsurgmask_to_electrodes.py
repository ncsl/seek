import argparse
import os
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from mne_bids import BIDSPath
from mne_bids.tsv_handler import _from_tsv
from nibabel.affines import apply_affine

sys.path.append("../../../")


def apply_mask_on_electrodes(electrodes_tsv, mask_img):
    inv_affine = npl.inv(
        mask_img.affine
    )  # Obtain inverse affine matrix to transform from xyz to voxel

    # convert electrode coordinates to voxels
    ch_names = electrodes_tsv["name"]
    elecmatrix = []
    for i in range(len(ch_names)):
        elecmatrix.append([electrodes_tsv[x][i] for x in ["x", "y", "z"]])
    elecmatrix = np.array(elecmatrix, dtype=float)

    # apply affine transform to put into Voxel space
    electrodes_vox = apply_affine(inv_affine, elecmatrix)

    # loop over channels and their voxel coordinates
    resected_contacts = dict()

    mask_indx = np.argwhere(mask_img.get_fdata())
    # print(np.where(mask_img.get_fdata() != 0))
    for ch_name, vox_coords in zip(ch_names, electrodes_vox):
        vox_coords = np.array(list(map(int, vox_coords)))
        # print(vox_coords)
        # print(mask_img.shape)
        # print(mask_img.get_fdata()[vox_coords])
        # if mask_img.get_fdata()[vox_coords[0], vox_coords[1], vox_coords[2]] != 0:
        #     resected_contacts.append(ch_name)

        for i in range(len(mask_indx)):
            if (np.abs(vox_coords[0] - mask_indx[i][0]) <= 3) & \
                    (np.abs(vox_coords[1] - mask_indx[i][1]) <= 5.0) & \
                    (np.abs(vox_coords[2] - mask_indx[i][2]) <= 3.1):
                resected_contacts[ch_name] = 1
    print(resected_contacts.keys())


if __name__ == "__main__":
    root = Path("/Users/adam2392/OneDrive - Johns Hopkins/epilepsy_bids/")
    sourcepath = root / "sourcedata"
    subject = "la02"
    space = "fs"
    acquisition = "seeg"
    session = "presurgery"

    # FreeSurfer file items
    fs_subj_dir = root / "derivatives" / "freesurfer" / subject
    fs_lut_fpath = sourcepath / "electrodes localized" / "FreeSurferColorLUT.txt"

    # default output electrodes tsv file path
    electrodes_tsv_fpath = BIDSPath(
        subject=subject,
        session=session,
        space=space,
        acquisition=acquisition,
        datatype="ieeg",
        root=root,
        suffix="electrodes",
        extension=".tsv",
    )

    mask_img_fpath = BIDSPath(
        subject=subject,
        session='postsurgery',
        datatype="anat",
        processing='slicer',
        space=space,
        suffix="mask",
        extension=".nii",
        root=root,
        check=False
    )

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-mri_xyzcoords_fpath",
        required=False,
        help="The output datafile with all the electrode points clustered.",
        default=sourcepath
                / "electrodes localized"
                / "stolk"
                / f"{subject}_elec_acpc_f.mat",
    )
    parser.add_argument(
        "-output_bids_electrodes_file",
        help="The output BIDS datafile for electrodes in tsv format.",
        required=False,
        default=electrodes_tsv_fpath.fpath,
    )
    parser.add_argument(
        "-fs_patient_dir",
        help="The freesurfer output directory.",
        required=False,
        default=fs_subj_dir,
    )
    parser.add_argument(
        "-fs_lut_fpath",
        help="The Freesurfer LUT.",
        required=False,
        default=fs_lut_fpath,
    )
    parser.add_argument("-mask_img_fpath", required=False, default=mask_img_fpath.fpath)
    args = parser.parse_args()

    # Extract arguments from parser
    mri_xyzcoords_fpath = args.mri_xyzcoords_fpath
    output_electrodes_tsv_fpath = args.output_bids_electrodes_file
    fs_patient_dir = args.fs_patient_dir
    fs_lut_fpath = args.fs_lut_fpath
    mask_img_fpath = args.mask_img_fpath

    # load in the T1 MRI image and its affine
    print(f'Loading in mask file {mask_img_fpath}')
    t1_img = nb.load(mask_img_fpath)

    # load in the electrode coordinates
    print(f'Loading file {output_electrodes_tsv_fpath}')
    electrodes_tsv = _from_tsv(output_electrodes_tsv_fpath)

    # extract subject id from bids sidecar electrodes fname
    subject_id = (
        os.path.basename(output_electrodes_tsv_fpath).split("_")[0].split("sub-")[1]
    )

    # apply atlas
    apply_mask_on_electrodes(electrodes_tsv, t1_img)
