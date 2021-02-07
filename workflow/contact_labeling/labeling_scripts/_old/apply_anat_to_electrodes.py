import argparse
import json
import os
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import numpy.linalg as npl
import scipy.io
from mne_bids import BIDSPath, get_entities_from_fname
from mne_bids.tsv_handler import _from_tsv, _to_tsv
from nibabel.affines import apply_affine

from seek_localize import read_dig_bids, label_elecs_anat, fs_lut_fpath

sys.path.append("../../../")
from seek.format.bids_conversion import _write_coordsystem_json
from seek.contacts.anat.label_chs_anat import (
    convert_fsmesh2mlab,
    label_elecs,
)
from seek.format.bids_conversion import (
    _update_electrodes_json,
    _update_electrodes_tsv,
)


def save_bids_files(coordsystem_fpath, mri_img_fpath):
    _write_coordsystem_json(fname=coordsystem_fpath, unit="mm", img_fname=mri_img_fpath)


def label_anatomy(root, elecs_fname, atlas_img_fname):
    # get paths
    entities = get_entities_from_fname(elecs_fname)
    elecs_path = BIDSPath(**entities, root=root, extension=".tsv")
    coordsystem_path = elecs_path.copy().update(suffix="coordsystem", extension=".json")

    # label anatomy
    elecs_df = label_elecs_anat(
        bids_path=elecs_path, img_fname=atlas_img_fname, fs_lut_fpath=fs_lut_fpath
    )

    # save to disc
    elecs_df.to_csv(elecs_path, sep="\t", index=None)

    # create sidecar electrodes json file
    electrodes_json_fpath = str(elecs_path).replace(".tsv", ".json")
    json_dict = {
        "destriuex": "Electrode annotation using Destriuex atlas with 196 brain regions.",
        "desikan-killiany": "Electrode annotation using DK atlas with 86 brain regions.",
    }
    electrodes_json = _update_electrodes_json(electrodes_json_fpath, **json_dict)

    return elecs_df, electrodes_json


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
    root = Path("/Users/adam2392/OneDrive - Johns Hopkins/epilepsy_bids/")
    sourcepath = root / "sourcedata"
    subject = "la02"
    space = "fs"
    acquisition = "seeg"
    session = "presurgery"

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
    # FreeSurfer file items
    fs_subj_dir = root / "derivatives" / "freesurfer" / subject
    fs_lut_fpath = sourcepath / "electrodes localized" / "FreeSurferColorLUT.txt"
    mri_img_fpath = BIDSPath(
        subject=subject,
        session=session,
        datatype="anat",
        space=space,
        suffix="T1w",
        extension=".nii",
        root=root,
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
    parser.add_argument("-mri_img_fpath", required=False, default=mri_img_fpath.fpath)
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
    if str(mri_xyzcoords_fpath).endswith(".tsv"):
        electrodes_tsv = _from_tsv(mri_xyzcoords_fpath)
    else:
        # read in mat file
        electrodes_tsv = read_label_coords(mri_xyzcoords_fpath)

    # extract subject id from bids sidecar electrodes fname
    subject_id = (
        os.path.basename(output_electrodes_tsv_fpath).split("_")[0].split("sub-")[1]
    )

    # bids_root
    freesurfer_dir = Path(fs_patient_dir).parent
    derivatives_dir = freesurfer_dir.parent
    bids_root = derivatives_dir.parent

    # save it to electrodes output tsv
    # write the output to a txt file
    with open(output_electrodes_tsv_fpath, "w", encoding="utf-8") as f:
        f.write("name\tx\ty\tz\n")
        for i, name in enumerate(electrodes_tsv.keys()):
            f.write(
                "%s\t%.6f\t%.6f\t%.6f\n"
                % (
                    name,
                    electrodes_tsv[name][0],
                    electrodes_tsv[name][1],
                    electrodes_tsv[name][2],
                )
            )

    # Output labeled .mat files with atlas, white matter, and brainmask information
    electrodes_tsv, electrodes_json = apply_atlas(
        bids_root, output_electrodes_tsv_fpath, inv_affine, fs_patient_dir, fs_lut_fpath
    )

    # save sidecar electrodes tsv
    _to_tsv(electrodes_tsv, output_electrodes_tsv_fpath)

    # save sidecar electrodes json
    with open(str(output_electrodes_tsv_fpath).replace(".tsv", ".json"), "w") as fout:
        json.dump(electrodes_json, fout, indent=4)

    # create a coordsystem JSON file
    output_coordsyste_fpath = str(output_electrodes_tsv_fpath).replace(
        "electrodes.tsv", "coordsystem.json"
    )
    _write_coordsystem_json(
        fname=output_coordsyste_fpath, unit="mm", img_fname=mri_img_fpath
    )
