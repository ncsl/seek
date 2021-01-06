"""Untested API for converting neuroimaging files to BIDS format."""
import json
import os
import platform
import tempfile
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import List, Union

import dicom2nifti
import nibabel as nb
import numpy as np
from mne.utils import run_subprocess
from mne_bids import write_anat
from mne_bids.tsv_handler import _from_tsv, _to_tsv
from mne_bids.write import (
    _write_json,
    _write_tsv,
)
from mne_bids.sidecar_updates import _update_sidecar
from mne_bids.path import get_entities_from_fname


def bids_validate(bids_root):
    """Run BIDS validator."""
    shell = False
    bids_validator_exe = ["bids-validator", "--config.error=41", "--conwfig.error=41"]
    if platform.system() == "Windows":
        shell = True
        exe = os.getenv("VALIDATOR_EXECUTABLE", "n/a")
        if "VALIDATOR_EXECUTABLE" != "n/a":
            bids_validator_exe = ["node", exe]

    def _validate(bids_root):
        cmd = bids_validator_exe + [bids_root]
        run_subprocess(cmd, shell=shell)

    return _validate(bids_root)


def _update_electrodes_tsv(electrodes_tsv_fpath, elec_labels_anat, atlas_depth):
    electrodes_tsv = _from_tsv(electrodes_tsv_fpath)

    if atlas_depth not in electrodes_tsv.keys():
        electrodes_tsv[atlas_depth] = ["n/a"] * len(elec_labels_anat)
    for i in range(len(elec_labels_anat)):
        ch_name = electrodes_tsv["name"][i]
        print(ch_name, elec_labels_anat[i])
        electrodes_tsv[atlas_depth][i] = elec_labels_anat[i]
    _to_tsv(electrodes_tsv, electrodes_tsv_fpath)

    return electrodes_tsv


def _update_electrodes_json(electrodes_json_fpath, **kwargs):
    electrodes_json_fpath = Path(electrodes_json_fpath)
    if not electrodes_json_fpath.exists():
        electrodes_json_fpath.parent.mkdir(parents=True, exist_ok=True)
        with open(electrodes_json_fpath, "w") as fout:
            sidecar_json = json.dump(kwargs, fout)
    else:
        for key, val in kwargs.items():
            _update_sidecar(electrodes_json_fpath, key, val)
        with open(electrodes_json_fpath, "r") as fin:
            sidecar_json = json.load(fin)
    return sidecar_json


def _write_coordsystem_json(
    fname: str,
    unit: str,
    img_fname: str = None,
    overwrite: bool = True,
    verbose: bool = True,
    coordsystem: str = None,
):
    system_description = (
        "FreeSurfer Coordinate System derived from the CT, or T1 MRI scan."
    )
    processing_description = "SEEK-algorithm (thresholding, cylindrical clustering and post-processing), or manual labeling of contacts using FieldTrip Toolbox."

    if img_fname is not None:
        # load in image and determine coordinate system
        img = nb.load(img_fname)
        axcodes = nb.orientations.aff2axcodes(img.affine)
        coordsystem_name = "".join(axcodes)
    else:
        warnings.warn(
            "Image filename not passed in... Defaulting to MRI coordinate system."
        )
        coordsystem_name = "MRI"

    if coordsystem is None:
        coordsystem = coordsystem_name

    fid_json = {
        "IntendedFor": os.path.basename(img_fname),
        "iEEGCoordinateSystem": coordsystem,  # MRI, Pixels, or ACPC
        "iEEGCoordinateUnits": unit,  # m (MNE), mm, cm , or pixels
        "iEEGCoordinateSystemDescription": system_description,
        "iEEGCoordinateProcessingDescription": processing_description,
        "iEEGCoordinateProcessingReference": "See DOI: https://zenodo.org/record/3542307#.XoYF9tNKhZI",
    }
    _write_json(fname, fid_json, overwrite, verbose)

    return fname


def _write_electrodes_tsv(
    fname: str,
    ch_names: Union[List, np.ndarray],
    coords: Union[List, np.ndarray],
    sizes: Union[List, np.ndarray] = None,
    overwrite: bool = False,
    verbose: bool = True,
):
    """
    Create an electrodes.tsv file and save it.

    Parameters
    ----------
    fname : str
        Filename to save the electrodes.tsv to.
    names :
    coords :
    sizes :
    overwrite : bool
        Defaults to False.
        Whether to overwrite the existing data in the file.
        If there is already data for the given `fname` and overwrite is False,
        an error will be raised.
    verbose :  bool
        Set verbose output to true or false.
    """
    if len(ch_names) != len(coords):
        raise RuntimeError(
            "Number of channel names should match "
            "number of coordinates passed in. "
            f"{len(ch_names)} names and {len(coords)} coords passed in."
        )

    x, y, z, names = list(), list(), list(), list()
    for name, coord in zip(ch_names, coords):
        x.append(coord[0])
        y.append(coord[1])
        z.append(coord[2])
        names.append(name)

    if sizes is None:
        sizes = ["n/a"] * len(ch_names)

    data = OrderedDict(
        [
            ("name", names),
            ("x", x),
            ("y", y),
            ("z", z),
            ("size", sizes),
        ]
    )

    print(f"Writin data to {fname}: ")
    print(data)

    _write_tsv(fname, data, overwrite=overwrite, verbose=verbose)
    return fname


def _convert_dicom_to_nifti(original_dicom_directory, output_fpath, verbose=True):
    if verbose:
        print(f"Reading in dicoms from {original_dicom_directory}")
        print(f"Saving to {output_fpath}")

    # try to run mrconvert and reorient to `LAS` direction
    output_dict = dicom2nifti.dicom_series_to_nifti(
        original_dicom_directory, output_fpath, reorient_nifti=False
    )
    nb_img = output_dict["NII"]
    image_input = output_dict["NII_FILE"]

    print("Orientation of nifti image: ", nb.aff2axcodes(nb_img.affine))
    img = nb.load(image_input)
    print("Reoriented image to ", nb.aff2axcodes(img.affine))

    return image_input


def convert_img_to_bids(image_input, bids_root, bids_fname, verbose=True):
    """Run Bids Conversion script to be updated.

    Performs BIDS conversion for Ct/T1/DTI/fMRI data.

    TODO: demo for DTI/FMRI
    """
    if verbose:
        print(f"bids_root is {bids_root}")
        print(f"Reading in image files from: {image_input}")

    # create temporary filepath to store the nifti file
    with tempfile.TemporaryDirectory() as tmpdir:
        output_fpath = Path(tmpdir, "tmp.nii").as_posix()

        if len([x for x in Path(image_input).glob("*.dcm")]) > 0:
            print("Converting dicom -> Nifti...")
            # try to run mrconvert and reorient to `LAS` direction
            # try:
            image_input = _convert_dicom_to_nifti(image_input, output_fpath)
            # except Exception as e:
            #     "mrconvert {params.CT_FOLDER} {output.CT_bids_fname};"
        else:
            print("Passed NIFTI image, so skipping mrconvert from dicom -> nifti...")
            image_input = str(image_input)

        print(image_input)
        # determine the BIDS identifiers
        params = get_entities_from_fname(bids_fname)
        subject = params["sub"]
        session = params["ses"]

        print("\n\nWriting now to BIDS...")
        # write to BIDS
        anat_dir = write_anat(
            bids_root,
            subject,
            t1w=image_input,
            session=session,
            overwrite=True,
            verbose=True,
        )

    # image_fpath = Path(anat_dir, bids_fname).as_posix()
    # img = nb.load(image_fpath)
    # print(img)
    # print("Reoriented image to ", nb.aff2axcodes(img.affine))
    # validate that this is bids valid dataset
    # try:
    #     bids_validate(bids_root)
    # except Exception as e:
    #     print(e)

    # Plot it
    # from nilearn.plotting import plot_anat
    # import matplotlib.pyplot as plt
    #
    # fig, axs = plt.subplots(3, 1)
    # for point_idx, label in enumerate(("LPA", "NAS", "RPA")):
    #     plot_anat(Path(anat_dir, bids_fname), axes=axs[point_idx], title=label)
    # plt.show()


def save_organized_elecdict_astsv(elecdict, output_fpath, size=None, img_fname=None):
    """Save organized electrode dict coordinates as a tsv file."""
    x, y, z, names = list(), list(), list(), list()
    coords = []
    for elec in elecdict.keys():
        for ch, ch_coord in elecdict[elec].items():
            x.append(ch_coord[0])
            y.append(ch_coord[1])
            z.append(ch_coord[2])
            names.append(ch)
            coords.append(ch_coord)
        if size is None:
            sizes = ["n/a"] * len(names)
        else:
            sizes = [size] * len(names)

    # write the tsv file
    _write_electrodes_tsv(output_fpath, names, coords, sizes)

    outputjson_fpath = output_fpath.replace("electrodes.tsv", "coordsystem.json")
    unit = "mm"
    # write accompanying coordinate system json file
    _write_coordsystem_json(outputjson_fpath, unit, img_fname)
