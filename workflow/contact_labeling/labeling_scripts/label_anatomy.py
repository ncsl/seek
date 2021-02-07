import argparse

from mne_bids import BIDSPath, get_entities_from_fname
from seek_localize import label_elecs_anat, fs_lut_fpath

from seek.format.bids_conversion import (
    _update_electrodes_json,
)


def label_anatomy(root, elecs_fname, atlas_img_fname):
    """Label anatomy script.

    Parameters
    ----------
    root : str
        Root of BIDS dataset.
    elecs_fname : str
        File path for electrodes file.
    atlas_img_fname : str
        File path for the atlas image volume.

    Returns
    -------
    elecs_df : pd.DataFrame
    electrodes_json : dict
    """
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


def main(root, elecs_fname, atlas_img_fname, output_elecs_fname):
    elecs_df, elecs_json = label_anatomy(root, elecs_fname, atlas_img_fname)

    # save to output specified location on disc
    elecs_df.to_csv(output_elecs_fname, sep="\t", index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-elecs_fpath",
        required=True,
        help="The input BIDS `*electrodes.tsv` file with all "
        "the electrode contacts labeled and localized. "
        "A corresponding `*coordsystem.json` file must "
        "also accompany it.",
    )
    parser.add_argument(
        "-img_fpath",
        required=True,
        help="The atlas image filepath that contains the annotated "
        "voxel points according to an atlas like Desikan-Killiany, "
        "or Destrieux."
        # help="The input BIDS image file (e.g. `*T1w.nii) corresponding "
        #      "to the image that contacts were localized in. "
        #      "Must be a Nifti file and must be BIDS-compliant. "
        #      "Ex name: `sub-01_ses-presurgery_space-fs_T1w.nii`."
    )
    parser.add_argument(
        "-output_elecs_fpath",
        help="The output BIDS datafile for electrodes in tsv format.",
        required=True,
    )
    parser.add_argument(
        "-bids_root",
        help="The root of the BIDS dataset.",
        required=False,
    )

    args = parser.parse_args()

    # Extract arguments from parser
    elecs_fpath = args.elecs_fpath
    img_fpath = args.img_fpath
    output_elecs_fpath = args.output_elecs_fpath
    bids_root = args.bids_root

    # run main script
    main(bids_root, elecs_fpath, img_fpath, output_elecs_fpath)
