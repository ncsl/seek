# -*- coding: utf-8 -*-

import json
import logging
import os
import re
import subprocess
import tempfile
from datetime import date, datetime
from functools import reduce
from typing import List, Optional

import nibabel
import nibabel as nb
import numpy as np
import pandas as pd

from seek.base.objects.baseneuroimage import Hemisphere
from seek.base.objects.baseneuroimage import RegionIndexMapping

"""
 -patient.apply_xfm
(xfm_dir='mri/transforms',
xfm_file='talairach.xfm',
source_file='elecs/TDT_elecs_all.mat',
target_file='elecs/
TDT_elecs_all_2tal.mat',
file_type='elecs')
- patient.apply_xfm
(xfm_dir='mri/transforms',
xfm_file='talairach.xfm',
source_file='Meshes/
rh_pial_trivert.mat',
target_file='Meshes/
rh_pial_trivert_2tal.mat',
file_type='surf')
- patient.apply_xfm
(xfm_dir='mri/transforms',
xfm_file='talairach.xfm',
source_file='Meshes/
lh_pial_trivert.mat',
target_file='Meshes/
lh_pial_trivert_2tal.mat',file_type='surf')"""


def apply_xfm_to_elecs(xfm_fpath, electrodes_fname, t1w_fpath):
    # First load the matrix in the xfm file:
    xfm_arr = np.genfromtxt(xfm_fpath, skip_header=5, delimiter=" ", comments=";")

    # get the vox2ras affine
    t1w = nb.load(t1w_fpath)
    vox2ras = t1w.affine
    vox2ras_tk = vox2ras

    # load in the electrodes data
    electrodes_tsv = pd.read_csv(electrodes_fname, delimiter="\t")
    ch_names = electrodes_tsv["name"].tolist()
    ch_coords = electrodes_tsv[["x", "y", "z"]].to_numpy(dtype=float)
    native_coords = []
    for name, coord in zip(ch_names, ch_coords):
        native_coord = _apply_xfm_to_native(xfm_arr, vox2ras, vox2ras_tk, coord)
        native_coords.append(native_coord)
    native_coords = np.array(native_coords)
    return native_coords


def _apply_xfm_to_native(xfm_arr, vox2ras, vox2ras_tk, input_vec):
    if len(input_vec) == 3:
        input_vec = np.concatenate((input_vec, [1]))

    # perform seriees of affine transformations to obtain native space coordinate
    input_vec_xfm = reduce(
        np.dot, [xfm_arr, vox2ras, np.linalg.inv(vox2ras_tk), input_vec.T]
    ).T
    return input_vec_xfm


def _contact_numbers_on_electrode(entry_ch, exit_ch):
    elec_name1, entry_ch_num = re.match("^([A-Za-z]+[']?)([0-9]+)$", entry_ch).groups()
    elec_name2, exit_ch_num = re.match("^([A-Za-z]+[']?)([0-9]+)$", exit_ch).groups()
    if elec_name1 != elec_name2:
        raise RuntimeError(
            "Electrode name for entry/exit point should be the same. "
            f"Your electrode names for entry/exit are: {entry_ch}, and {exit_ch}."
        )

    # create integer list
    elec_contact_nums = np.arange(int(entry_ch_num), int(exit_ch_num) + 1, dtype=int)
    return elec_contact_nums


def group_contacts(elec_in_brain):
    """
    Group individual contacts by the electrode to which they correspond.

    Sorts the contacts using the corresponding labels.

    Parameters
    ----------
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
            key=lambda x: int(re.findall(r"\d+", x[0])[0]),
        )
        labeled_contacts[elec] = dict(sorted_chans)

    return labeled_contacts


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):  # This is the fix
            return obj.tolist()
        elif isinstance(obj, (datetime, date)):
            return obj.isoformat()
        return json.JSONEncoder.default(self, obj)


def ensure_string(arg):
    if not (isinstance(arg, str)):
        if arg is None:
            return ""
        else:
            return ensure_list(arg)[0]
    else:
        return arg


def ensure_list(arg):
    if not (isinstance(arg, list)):
        try:  # if iterable
            if isinstance(arg, (str, dict)):
                arg = [arg]
            else:
                arg = list(arg)
        except BaseException:  # if not iterable
            arg = [arg]
    return arg


def generate_region_labels(n_regions, labels=[], str=". ", numbering=True):
    if len(labels) == n_regions:
        if numbering:
            return np.array(
                [
                    str.join(["%d", "%s"]) % tuple(l)
                    for l in zip(range(n_regions), labels)
                ]
            )
        else:
            return labels
    else:
        return np.array(["%d" % l for l in range(n_regions)])


def pial_to_verts_and_triangs(pial_surf) -> (np.ndarray, np.ndarray):
    """
    Convert pial surface file to vertices and triangles

    Parameters
    ----------
    pial_surf :

    Returns
    -------

    """
    tmpdir = tempfile.TemporaryDirectory()
    pial_asc = os.path.join(tmpdir.name, os.path.basename(pial_surf + ".asc"))
    subprocess.run(["mris_convert", pial_surf, pial_asc])

    with open(pial_asc, "r") as f:
        f.readline()
        nverts, ntriangs = [int(n) for n in f.readline().strip().split(" ")]

    vertices = np.genfromtxt(
        pial_asc, dtype=float, skip_header=2, skip_footer=ntriangs, usecols=(0, 1, 2)
    )
    triangles = np.genfromtxt(
        pial_asc, dtype=int, skip_header=2 + nverts, usecols=(0, 1, 2)
    )
    assert vertices.shape == (nverts, 3)
    assert triangles.shape == (ntriangs, 3)

    completed_process = subprocess.run(
        ["mris_info", pial_surf], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    mris_info = completed_process.stdout.decode("ascii")
    c_ras_list = extract_vector(mris_info, "c_(ras)")
    assert c_ras_list is not None
    vertices[:, 0:3] += np.array(c_ras_list)

    return vertices, triangles


def read_cortical_region_mapping(
    label_direc: os.PathLike, hemisphere: Hemisphere, fs_to_conn: RegionIndexMapping
) -> np.ndarray:
    """
    Reads the cortical region mapping file.

    :param label_direc: Where the annotation label directory is
    :param hemisphere: (Hemisphere) enumerator
    :param fs_to_conn: (RegionIndexMapping)
    :return:
    """
    filename = os.path.join(label_direc, hemisphere.value + ".aparc.annot")
    region_mapping, _, _ = nibabel.freesurfer.io.read_annot(filename)
    region_mapping = region_mapping - 1
    region_mapping[region_mapping == -2] = 0  # Unknown regions in hemispheres

    # $FREESURFER_HOME/FreeSurferColorLUT.txt describes the shift
    if hemisphere == Hemisphere.lh:
        region_mapping += FS_LUT_LH_SHIFT
    else:
        region_mapping += FS_LUT_RH_SHIFT

    fs_to_conn_fun = np.vectorize(lambda n: fs_to_conn.source_to_target(n))
    region_mapping = fs_to_conn_fun(region_mapping)

    return region_mapping


def extract_vector(string: str, name: str) -> Optional[List[float]]:
    r"""
    Extract numerical vector from a block of text. The vector has to be on a single line with the format:
    <name> : (x0, x1, x2 [,...])
    If the vector in the correct format is missing, return None.

    >>> extract_vector("EXAMPLE\na: (1.0, 2.0, 3.0)\nb: (0.0, 0.0, 0.0)", "a")
    [1.0, 2.0, 3.0]

    >>> extract_vector("EMPTY", "a") is None
    True
    """

    for line in string.split("\n"):
        match = re.match(
            r"""^\s*
                         (.+?)              # name
                         \s*:\s*            # separator
                         \(([0-9.,\s-]+)\)   # vector: (x0, x1, ....)
                         \s*$""",
            line,
            re.X,
        )

        if match and match.group(1) == name:
            try:
                vector = [float(x) for x in match.group(2).split(",")]
                return vector
            except ValueError:
                pass

    return None


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
