"""Untested API for converting neuroimaging files to BIDS format."""
import os
import platform
import tempfile
from pathlib import Path

import dicom2nifti
import nibabel as nb
from mne.utils import run_subprocess
from mne_bids import write_anat
from mne_bids.utils import _parse_bids_filename


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


def _convert_dicom_to_nifti(original_dicom_directory, output_fpath, verbose=True):
    if verbose:
        print(f"Reading in dicoms from {original_dicom_directory}")
        print(f"Saving to {output_fpath}")

    # try to run mrconvert and reorient to `LAS` direction
    output_dict = dicom2nifti.dicom_series_to_nifti(
        original_dicom_directory, output_fpath, reorient_nifti=True
    )
    nb_img = output_dict["NII"]
    image_input = output_dict["NII_FILE"]

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
        params = _parse_bids_filename(bids_fname, verbose=True)
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


if __name__ == "__main__":
    # bids root to write BIDS data to
    bids_root = Path("/home/adam2392/hdd2/data/")
    # bids_root = Path("/home/adam2392/Documents/eztrack/data/bids_layout/")
    # bids_root = Path("/Users/adam2392/Documents/eztrack/data/bids_layout/")

    # path to original source data
    center = "clevelandnl"
    source_path = Path(bids_root / "sourcedata" / center)

    # HACK: get all subject ids within sourcedata
    subject_ids = natsorted(
        [
            x.name
            for x in source_path.iterdir()
            if not x.as_posix().startswith(".")
            if x.is_dir()
        ]
    )[0:1]

    # define BIDS identifiers
    task = "monitor"
    session = "seizure"
