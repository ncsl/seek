# -*- coding: utf-8 -*-

import logging
import os

import nibabel as nb
import numpy as np


def label_volume_centers(label_volume, output_tsv):
    """
    A function to take in freesurfer output aparc+aseg.mgz 
    which is a volume file to label each of the volume centers.

    To examine the Desikan-Killiany atlas, look at: aparc+aseg.mgz
    - for Destrieux, look at: aparc.a2009s+aseg.mgz
    """
    log = logging.getLogger("label_volume_centers")
    log.info("reading %r", label_volume)

    vol = nb.load(label_volume)

    log.info("computing centers")
    centers = list(compute_label_volume_centers(vol))

    log.info("loading FS LUT")
    lut_path = os.path.join(os.environ["FREESURFER_HOME"], "FreeSurferColorLUT.txt")
    lut_map = build_fs_label_name_map(lut_path)

    log.info("writing result to %r", output_tsv)
    with open(output_tsv, "w") as fd:
        for val, (x, y, z) in centers:
            val_ = lut_map[val] if lut_map else val
            fd.write("%f\t%f\t%f\t%s\n" % (x, y, z, val_))


def compute_label_volume_centers(label_volume, affine=None):
    """
    Get the centers for each unique value within the volume.

    For FreeSurfer this can range from any of the values they specify in their
    freesurfer directory's LookUp Tables (LUT) (e.g. FreeSurferColorLUT.txt). This
    allows you to map the numerical segmentation, parcellation to an anatomical region.

    """

    try:
        vol = label_volume.get_data()
        aff = label_volume.affine
    except:
        vol = label_volume
        aff = affine
    for val in np.unique(vol):
        xyz = vol_val_xyz(vol, aff, val)
        x, y, z = xyz.mean(axis=0)
        yield val, (x, y, z)


def vol_val_xyz(vol, aff, val):
    """
    Get the xyz position(s) in the volume, for a specific volume value.
    """
    vox_idx = np.argwhere(vol == val)
    xyz = aff.dot(np.c_[vox_idx, np.ones(vox_idx.shape[0])].T)[:3].T
    return xyz


def build_fs_label_name_map(lut_path):
    """
    Build's freesurfer label name mapping from values, to names and their
    corresponding lut values.
    """
    lut = {}
    with open(lut_path, "r") as fd:
        # read in path line by line
        for line in fd.readlines():
            # skip comment lines
            if not line[0] == "#" and line.strip():
                val, name, _, _, _, _ = line.strip().split()
                lut[int(val)] = name
    return lut


def transform(coords, src_img, dest_img, transform_mat):
    import subprocess

    coords_str = " ".join([str(x) for x in coords])

    cp = subprocess.run(
        "echo %s | img2imgcoord -mm -src %s -dest %s -xfm %s"
        % (coords_str, src_img, dest_img, transform_mat),
        shell=True,
        stdout=subprocess.PIPE,
    )
    transformed_coords_str = cp.stdout.decode("ascii").strip().split("\n")[-1]
    return np.array([float(x) for x in transformed_coords_str.split(" ") if x])
