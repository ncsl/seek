import argparse
import json
import os
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import numpy.linalg as npl
import scipy.io
from mne_bids.tsv_handler import _from_tsv, _to_tsv
from nibabel.affines import apply_affine

sys.path.append("../../../")

from seek.contacts.anat.label_chs_anat import (
    convert_fsmesh2mlab,
    label_elecs,
)
from seek.format.bids_conversion import (
    _update_electrodes_json,
    _update_electrodes_tsv,
)


def apply_wm_and_brainmask(final_centroids_xyz, atlasfilepath, wmpath, bmpath):
    """
    Apply white matter and brainmask labels to final centroid output and save
    in .mat files.

    Parameters
    ----------
        final_centroids_xyz: dict(str: dict(str: ndarray))
            Dictionary of predicted centroids in xyz (mm) coordinates.

        atlasfilepath: str
            Path to .txt file to save the xyz coordinates of centroids.

        wmpath: str
            Path to white matter mask file.

        bmpath: str
            Path to brain matter mask file.

    Returns
    -------
        anatomy: ndarray
            Anatomy matrix with columns of coordinates, anatomical label, and
            channel label.
    """
    dat = scipy.io.loadmat(atlasfilepath)
    elecmatrix = dat["elecmatrix"]
    anatomy_orig = dat["anatomy"]
    eleclabels = dat["eleclabels"]

    # Load white matter and brain masks
    wm_img = nb.load(wmpath)
    wm_dat = wm_img.get_data()
    bm_img = nb.load(bmpath)
    bm_dat = bm_img.get_data()

    affine = npl.inv(bm_img.affine)

    wm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    bm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    # Add two columns in anatomy to store boolean for whether voxel is
    # white matter or brain matter
    anatomy = np.zeros((anatomy_orig.shape[0], anatomy_orig.shape[1] + 2), dtype=object)
    for i, label in enumerate(anatomy_orig):
        chan = str(label[0][0]).strip()
        for elec in final_centroids_xyz:
            if chan in final_centroids_xyz[elec].keys():
                pt = apply_affine(affine, final_centroids_xyz[elec][chan])
                wm_label[i] = wm_dat[list(map(int, pt))] > 0
                bm_label[i] = bm_dat[list(map(int, pt))] > 0
    anatomy[:, : anatomy_orig.shape[1]] = anatomy_orig
    anatomy[:, anatomy_orig.shape[1]] = wm_label
    anatomy[:, anatomy_orig.shape[1] + 1] = bm_label

    save_dict = {"elecmatrix": elecmatrix, "anatomy": anatomy, "eleclabels": eleclabels}
    scipy.io.savemat(atlasfilepath, mdict=save_dict)
    return anatomy


def apply_atlas(bids_root, electrodes_tsv_fpath, inv_affine, fspatdir, fs_lut_fpath):
    """
    Map centroids to an atlas (e.g. Desikan-Killiany, Destriuex) and apply
    white matter and brain masks to label centroids as white matter or out of
    the brain.

    Parameters
    –---------
        fspatdir: str
            Path to freesurfer directory.

        fs_lut_fpath: str
            Path to fs_lut file.

    Returns
    -------
        elec_labels_destriuex: dict(str: ndarray)
            array of contacts labeled with Destriuex atlas.

        elec_labels_DKT: dict(str: ndarray)
            array of contacts labeled with Desikan-Killiany atlas.
    """
    patid = os.path.basename(os.path.normpath(fspatdir))

    # Apply Atlases, white matter mask, and brainmask
    convert_fsmesh2mlab(subj_dir=os.path.abspath(os.path.dirname(fspatdir)), subj=patid)

    # load electrodes tsv
    electrodes_tsv = _from_tsv(electrodes_tsv_fpath)
    ch_names = electrodes_tsv["name"]
    elecmatrix = []
    for i in range(len(ch_names)):
        elecmatrix.append([electrodes_tsv[x][i] for x in ["x", "y", "z"]])
    elecmatrix = np.array(elecmatrix, dtype=float)

    # apply affine transform to put into Voxel space
    elecmatrix = apply_affine(inv_affine, elecmatrix)

    # anatomically label
    elec_labels_anat_destriuex = label_elecs(
        bids_root,
        ch_names,
        elecmatrix,
        subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
        subj=patid,
        hem="lh",
        fs_lut_fpath=fs_lut_fpath,
        atlas_depth="destriuex",
    )

    elec_labels_anat_dk = label_elecs(
        bids_root,
        ch_names,
        elecmatrix,
        subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
        subj=patid,
        hem="lh",
        fs_lut_fpath=fs_lut_fpath,
        atlas_depth="desikan-killiany",
    )

    # add atlas labeling to electrodes tsv data
    atlas_depth = "destriuex"
    electrodes_tsv = _update_electrodes_tsv(
        electrodes_tsv_fpath, elec_labels_anat_destriuex, atlas_depth
    )
    atlas_depth = "desikan-killiany"
    electrodes_tsv = _update_electrodes_tsv(
        electrodes_tsv_fpath, elec_labels_anat_dk, atlas_depth
    )

    # create sidecar electrodes json file
    electrodes_json_fpath = str(electrodes_tsv_fpath).replace(".tsv", ".json")
    json_dict = {
        "destriuex": "Electrode annotation using Destriuex atlas with 196 brain regions.",
        "desikan-killiany": "Electrode annotation using DK atlas with 86 brain regions.",
    }
    electrodes_json = _update_electrodes_json(electrodes_json_fpath, **json_dict)

    return electrodes_tsv, electrodes_json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mri_xyzcoords_fpath",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument(
        "output_bids_electrodes_file",
        help="The output BIDS datafile for electrodes in tsv format.",
    )
    parser.add_argument("fs_patient_dir", help="The freesurfer output diretroy.")
    parser.add_argument("fs_lut_fpath", help="The Freesurfer LUT.")
    parser.add_argument("mri_img_fpath", default=None)
    args = parser.parse_args()

    # Extract arguments from parser
    mri_xyzcoords_fpath = args.mri_xyzcoords_fpath
    output_electrodes_tsv_fpath = args.output_bids_electrodes_file
    fs_patient_dir = args.fs_patient_dir
    fs_lut_fpath = args.fs_lut_fpath
    mri_img_fpath = args.mri_img_fpath

    # load in the T1 MRI image and its affine
    t1_img = nb.load(mri_img_fpath)
    inv_affine = npl.inv(
        t1_img.affine
    )  # Obtain inverse affine matrix to transform from xyz to CT voxel

    # load in the electrode coordinates
    electrodes_tsv = _from_tsv(mri_xyzcoords_fpath)

    # extract subject id from bids sidecar electrodes fname
    subject_id = (
        os.path.basename(output_electrodes_tsv_fpath).split("_")[0].split("sub-")[1]
    )

    # bids_root
    freesurfer_dir = Path(fs_patient_dir).parent
    derivatives_dir = freesurfer_dir.parent
    bids_root = derivatives_dir.parent

    # Output labeled .mat files with atlas, white matter, and brainmask information
    electrodes_tsv, electrodes_json = apply_atlas(
        bids_root, mri_xyzcoords_fpath, inv_affine, fs_patient_dir, fs_lut_fpath
    )

    # save sidecar electrodes tsv
    _to_tsv(electrodes_tsv, output_electrodes_tsv_fpath)

    # save sidecar electrodes json
    with open(output_electrodes_tsv_fpath.replace(".tsv", ".json"), "w") as fout:
        json.dump(electrodes_json, fout)
