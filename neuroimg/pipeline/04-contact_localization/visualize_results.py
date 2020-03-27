import argparse
import math
import os
import sys

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import numpy.linalg as npl
import seaborn as sns
from mne_bids.tsv_handler import _from_tsv
from sklearn.decomposition import PCA

sys.path.append("../../../")

from neuroimg.base.utils.utils import MatReader, group_contacts


def summary_PCA_plots(figurefilepath, final_xyz_centroids, validation_xyz_centroids):
    """
    Construct PCA plots to visually show deviations in predicted centroids
    and manually labeled centroids.

    Parameters
    ----------
        figurefilepath: str
            File path to save the image to

        final_xyz_centroids: dict()
            Dictionary of the centroids for each contact along an electrode
            found by clustering algorithm.

        validation_xyz_centroids: dict()
            Dictionary of manually labeled centroids
    """
    if not os.path.exists(os.path.dirname(figurefilepath)):
        os.mkdir(os.path.dirname(figurefilepath))

    # convert these to electrodes dictionary
    final_xyz_centroids_elecs = group_contacts(final_xyz_centroids)
    validation_xyz_centroids_elecs = group_contacts(validation_xyz_centroids)

    # apply analysis to electrodes
    numelectrodes = len(final_xyz_centroids_elecs.keys())
    print(f"Plotting results for {numelectrodes} electrodes.")

    sns.set_context("paper", font_scale=1.1)
    fig, axs = plt.subplots(
        nrows=int(np.ceil(numelectrodes / 2)), ncols=2, figsize=(20, 20), dpi=100
    )
    axs = axs.flatten()

    for i, elec in enumerate(final_xyz_centroids_elecs.keys()):
        if elec not in validation_xyz_centroids_elecs.keys():
            print(f"Electrode {elec} not in validation centroids.")
            continue

        final_elec_coord_arr = final_xyz_centroids_elecs[elec].values()
        validation_elec_coord_arr = validation_xyz_centroids_elecs[elec].values()
        pred_centroids = np.array(list(final_elec_coord_arr))
        val_centroids = np.array(list(validation_elec_coord_arr))

        if pred_centroids.ndim == 1:
            pred_centroids = pred_centroids[:, np.newaxis]
        if val_centroids.ndim == 1:
            val_centroids = val_centroids[:, np.newaxis]

        # perform PCA onto a line
        pca = PCA(n_components=1).fit(val_centroids)
        pred_pca = pca.transform(pred_centroids)
        val_pca = pca.transform(val_centroids)

        axs[i].scatter(val_pca[:, 0], np.zeros_like(val_pca[:, 0]), label="expected")
        axs[i].scatter(
            pred_pca[:, 0],
            np.zeros_like(pred_pca[:, 0]),
            label="observed",
            marker="x",
            c="r",
        )
        axs[i].set(
            title=f"PCA Validation of Centroids (Electrode: {elec})",
            xlabel=f"PC 1st Coordinates in Voxels along Electrode {elec}",
            ylabel=f"PC 2nd Coordinates in Voxels along Electrode {elec}",
            ylim=[-0.005, 0.005],
        )
        axs[i].legend()
        for j, chan in enumerate(final_xyz_centroids_elecs[elec]):
            axs[i].annotate(chan, (pred_pca[j, 0], 0.0005), size=9)

        # set plot paramsx
    fig.tight_layout()
    plt.savefig(figurefilepath, box_inches="tight")
    plt.close(fig)

    return fig, axs


def l2_error(figurefilepath, final_centroids, validation_centroids):
    """
    Function that computes the Euclidean distance (L2) between centroids
    computed by algorithm and the validation_xyz_centroids data, which are centroids
    manually inputted by user.

    Parameters
    ----------
    final_centroids: dict()
        dictionary of properly labeled centroids grouped by electrode
        in CT voxels.

    validation_centroids: dict()
        electrode coordinates in CT voxels

    Returns
    -------
        Error in each channel (accurate prediction is between 0-3).
    """
    numchans = len(final_centroids.keys())

    errors_per_channel = {}
    for ch_name in final_centroids.keys():
        if ch_name not in validation_centroids.keys():
            print(f"Channel {ch_name} not in validation centroids.")
            continue

        val = validation_centroids[ch_name]  # Coordinates detected by algorithm
        pred = final_centroids[ch_name]
        abs_error = npl.norm(pred - val)
        errors_per_channel[ch_name] = abs_error

    # convert to errors per electrode
    errors_per_electrode = group_contacts(errors_per_channel)
    elec_names = errors_per_electrode.keys()
    numelecs = len(elec_names)

    sns.set(font_scale=1.1)
    fig, axs = plt.subplots(
        nrows=int(np.ceil(numelecs / 2)), ncols=2, figsize=(15, 15), dpi=100
    )
    axs = axs.flatten()

    ymin, ymax = 0, 20
    for i, elec_name in enumerate(errors_per_electrode.keys()):
        ch_names = errors_per_electrode[elec_name].keys()

        y_pos = np.arange(len(errors_per_electrode[elec_name]))
        axs[i].bar(
            y_pos,
            np.array(list(errors_per_electrode[elec_name].values())),
            align="center",
            alpha=0.9,
        )
        axs[i].set(
            title=f"Absolute Error By Channel in Electrode {elec_name}",
            xlabel="Channel",
            xticks=y_pos,
            xticklabels=list(ch_names),
            ylabel="Euclidean Distance (mm)",
            ylim=[ymin, ymax],
        )
        axs[i].set_xlabel("Channel")
        axs[i].set_ylabel("Euclidean Distance (mm)")
        axs[i].set_ylim([ymin, ymax])
    fig.tight_layout()

    plt.savefig(figurefilepath, box_inches="tight")

    return errors_per_channel, fig, axs


def load_data(elecfile):
    """
    Load electrodes file.

    Parameters
    ----------
    elecfile: str
        Path to space-delimited text file of contact labels and contact
        coordinates in mm space.

    Returns
    -------
        ct_img: NiBabel image object
            NiBabel image object of CT scan input.

        brainmask_ct: NiBabel image object
            NiBabel image object of brain mask in CT.

        elecinitfile: dict()
            A dictionary of contact coordinates in mm space. Keys are
            individual contact labels, and values are the corresponding
            coordinates in mm space.
    """

    elec_coords_mm = {}

    if elecfile.endswith(".txt"):
        with open(elecfile) as f:
            for l in f:
                row = l.split()
                elec_coords_mm[row[0]] = np.array(list(map(float, row[1:])))
    elif elecfile.endswith(".tsv"):
        electrodes_tsv = _from_tsv(elecfile)
        ch_names = electrodes_tsv["name"]
        elecmatrix = []
        for i in range(len(ch_names)):
            elecmatrix.append([electrodes_tsv[x][i] for x in ["x", "y", "z"]])
        elecmatrix = np.array(elecmatrix, dtype=float)

        for ch_name, coord in zip(ch_names, elecmatrix):
            elec_coords_mm[ch_name] = coord
    else:
        matreader = MatReader()
        data = matreader.loadmat(elecfile)
        elec_coords_mm = data["data"]

    print(f"\n\nFinished loading electrode locations data ({len(elec_coords_mm)}: ")
    print(list(elec_coords_mm.keys()))
    return elec_coords_mm


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "final_centroids_xyz_fpath", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "manual_centroids_xyz_fpath", help="The CT image volume in its original space."
    )
    parser.add_argument("ctfile", help="The CT img file")
    parser.add_argument("fsdir", help="The freesurfer output diretroy.")
    parser.add_argument("summary_pcaplots_filepath")
    parser.add_argument("summary_l2errorplots_filepath")
    args = parser.parse_args()

    # extract arguments from parser
    final_centroids_xyz_fpath = args.final_centroids_xyz_fpath
    manual_centroids_xyz_fpath = args.manual_centroids_xyz_fpath
    ctfile = args.ctfile
    fsdir = args.fsdir
    pcafig_filepath = args.summary_pcaplots_filepath
    l2fig_filepath = args.summary_l2errorplots_filepath

    # Load image data and get the image volumes
    ct_img = nb.load(ctfile)

    # load elecs data
    print("LOADING ELECS DATA FROM: ", final_centroids_xyz_fpath)
    elec_finalcoords_mm = load_data(final_centroids_xyz_fpath)

    print("LOADING ORIG ELECS DATA FROM: ", manual_centroids_xyz_fpath)
    elec_coords_mm = load_data(manual_centroids_xyz_fpath)

    """   SUMMARY PLOTS  """
    # create generated figures to check
    fig_pca, axs_pca = summary_PCA_plots(
        pcafig_filepath, elec_finalcoords_mm, elec_coords_mm
    )
    l2_errors, fig_l2, axs_l2 = l2_error(
        l2fig_filepath, elec_finalcoords_mm, elec_coords_mm
    )

    """   EUCLIDEAN DISTANCE ERROR CHECK """
    l2_errorstxt_path = os.path.join(
        os.path.dirname(pcafig_filepath), "euclidean_mm_distance_errors.txt"
    )
    with open(l2_errorstxt_path, "w") as f:
        for chan in l2_errors:
            err = l2_errors[chan]
            f.write("%s %.6f\n" % (chan, err))

    # Sum over all individual errors
    total_err = 0
    for chan in l2_errors:
        total_err += l2_errors[chan]

    print(f"Total L2 error in analysis: {total_err}")
