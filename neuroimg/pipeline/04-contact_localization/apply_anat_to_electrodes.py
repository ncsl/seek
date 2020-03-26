import argparse
import os
import sys

import nibabel as nb
import numpy as np
import numpy.linalg as npl
import scipy.io
from nibabel.affines import apply_affine
from mne_bids.tsv_handler import _from_tsv
from neuroimg.base.utils.io import load_elecs_data

sys.path.append("../../../")

from neuroimg.localize_contacts.freecog_labeling.utils import (
    convert_fsmesh2mlab,
    label_elecs,
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


def apply_atlas(fspatdir, destrieuxfilepath, dktfilepath, fs_lut_fpath):
    """
    Map centroids to an atlas (e.g. Desikan-Killiany, Destriuex) and apply
    white matter and brain masks to label centroids as white matter or out of
    the brain.

    Parameters
    –---------
        fspatdir: str
            Path to freesurfer directory.

        destrieuxfilepath: str
            Path to destrieux atlas for patient.

        dktfilepath: str
            Path to Desikan-Killiany atlas for patient.

        fs_lut_fpath: str
            Path to fs_lut file.

    Returns
    -------
        elec_labels_destriuex: dict(str: ndarray)
            array of contacts labeled with Destriuex atlas.

        elec_labels_DKT: dict(str: ndarray)
            array of contacts labeled with Desikan-Killiany atlas.
    """
    destriuexname = os.path.splitext(os.path.basename(destrieuxfilepath))[0]
    dktname = os.path.splitext(os.path.basename(dktfilepath))[0]

    patid = os.path.basename(os.path.normpath(fspatdir))

    # Apply Atlases, white matter mask, and brainmask
    convert_fsmesh2mlab(subj_dir=os.path.abspath(os.path.dirname(fspatdir)), subj=patid)
    elec_labels_destriuex = label_elecs(
        subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
        subj=patid,
        hem="lh",
        fs_lut_fpath=fs_lut_fpath,
        elecfile_prefix=destriuexname,
        atlas_depth="destriuex",
    )
    elec_labels_DKT = label_elecs(
        subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
        subj=patid,
        hem="lh",
        fs_lut_fpath=fs_lut_fpath,
        elecfile_prefix=dktname,
        atlas_depth="desikan-killiany",
    )
    return elec_labels_destriuex, elec_labels_DKT


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "clustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument("clustered_voxels_file", help="the voxels output datafile")
    parser.add_argument(
        "bids_electrodes_file",
        help="The output BIDS datafile for electrodes in tsv format.",
    )
    parser.add_argument("fs_patient_dir", help="The freesurfer output diretroy.")
    parser.add_argument("fs_lut_fpath", help="The Freesurfer LUT.")
    parser.add_argument("wm_native_file", default=None)
    args = parser.parse_args()

    # Extract arguments from parser
    clustered_points_file = args.clustered_points_file
    clustered_voxels_file = args.clustered_voxels_file
    electrodes_tsv_fpath = args.bids_electrodes_file
    fs_patient_dir = args.fs_patient_dir
    fs_lut_fpath = args.fs_lut_fpath

    # load in the electrode coordinates
    electrodes_tsv = _from_tsv(clustered_points_file)
    voxels_electrodes_dict = load_elecs_data(clustered_voxels_file)

    # Output labeled .mat files with atlas, white matter, and brainmask information
    elec_labels_destriuex, elec_labels_DKT = apply_atlas(
        fs_patient_dir, destrieuxfilepath, dktfilepath, fs_lut_fpath
    )

    """ SAVE CLUSTERED VOXELS AND POINTS AS TXT FILES WITH CHANNELS PER ROW """
    scipy.io.savemat(orgclustered_voxels_file, mdict={"data": final_centroids_voxels})
    scipy.io.savemat(orgclustered_points_file, mdict={"data": final_centroids_xyz})
