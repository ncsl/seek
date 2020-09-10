import argparse
import subprocess
import sys
from collections import OrderedDict

import numpy as np
from mne_bids.tsv_handler import _to_tsv, _from_tsv

sys.path.append("../../../")

from seek.utils.io import load_elecs_data
from seek.format.bids_conversion import _write_coordsystem_json


def transform(coords, src_img, dest_img, transform_mat, coordinate_type="mm"):
    coords_str = " ".join([str(x) for x in coords])

    cp = subprocess.run(
        f"echo %s | img2imgcoord -{coordinate_type} -src %s -dest %s -xfm %s"
        % (coords_str, src_img, dest_img, transform_mat),
        shell=True,
        stdout=subprocess.PIPE,
    )

    transformed_coords = cp.stdout.decode("ascii").strip().split("\n")[-1]
    # print(transformed_coords)
    return np.array([float(x) for x in transformed_coords.split(" ") if x])


def read_label_coords(elecfile):
    labels = []
    labelsxyz = []

    print("Reading ", elecfile)

    with open(elecfile, "r") as f:
        for _row in f:
            row = _row.split(" ")
            labels.append(row[0])
            labelsxyz.append([float(x) for x in row[1:]])

    return labels, labelsxyz


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ct_nifti_img", help="The CT image volume in its original space."
    )
    parser.add_argument("mri_nifti_img", help="Brain MRI image space.")
    parser.add_argument(
        "mapping_transformation_file", help="The mapping transformation file."
    )
    parser.add_argument(
        "ct_xyzcoords_fpath",
        help="The input datafile with all the electrode points in their coordinate space.",
    )
    parser.add_argument(
        "mri_xyzcoords_fpath",
        help="The output datafile for electrodes mapped to correct coords.",
    )
    parser.add_argument("--coordinate_type", default="mm")
    args = parser.parse_args()

    # extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    mri_nifti_img = args.mri_nifti_img
    mapping_transformation_file = args.mapping_transformation_file
    ct_coords_fpath = args.ct_xyzcoords_fpath
    mri_coords_fpath = args.mri_xyzcoords_fpath
    coordinate_type = args.coordinate_type

    # read in electrodes file
    if ct_coords_fpath.endswith(".tsv"):
        electrodes_tsv = _from_tsv(ct_coords_fpath)
    elif ct_coords_fpath.endswith(".txt"):
        electrodes_tsv = load_elecs_data(ct_coords_fpath)
        names = []
        x = []
        y = []
        z = []
        for ch_name, coord in electrodes_tsv.items():
            names.append(ch_name)
            x.append(coord[0])
            y.append(coord[1])
            z.append(coord[2])
        sizes = ["n/a"] * len(z)
        electrodes_tsv = OrderedDict(
            [
                ("name", names),
                ("x", x),
                ("y", y),
                ("z", z),
                ("size", sizes),
            ]
        )
    else:
        raise RuntimeError("CT Coords filepath needs to be txt or tsv.")

    # construct array of contact names and xyz coordinates
    labels = np.array(electrodes_tsv["name"])
    labelsxyz = []
    for i in range(len(labels)):
        labelsxyz.append([electrodes_tsv[x][i] for x in ["x", "y", "z"]])
    labelsxyz = np.array(labelsxyz)
    assert labelsxyz.shape[0] == len(labels)

    # run coordinate transformation
    modified_coords = np.array(
        [
            transform(
                coords,
                ct_nifti_img,
                mri_nifti_img,
                mapping_transformation_file,
                coordinate_type=coordinate_type,
            )
            for coords in labelsxyz
        ]
    )

    # resave the coordinate transformed data
    for i, (x, y, z) in enumerate(modified_coords):
        electrodes_tsv["x"][i] = x
        electrodes_tsv["y"][i] = y
        electrodes_tsv["z"][i] = z
    _to_tsv(electrodes_tsv, mri_coords_fpath)

    # resave the coordinate system file
    mri_coordsystem_fpath = mri_coords_fpath.replace(
        "electrodes.tsv", "coordsystem.json"
    )
    unit = "mm"
    img_fname = mri_nifti_img
    _write_coordsystem_json(mri_coordsystem_fpath, unit, img_fname=img_fname)
