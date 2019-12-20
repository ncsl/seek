import argparse
import math
import os
import sys

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import numpy.linalg as npl
import seaborn as sns
from sklearn.decomposition import PCA

sys.path.append("../../../")

from neuroimg.base.utils.utils import MatReader
from neuroimg.localize_contacts.electrode_clustering.mask import MaskVolume


def summary_PCA_plots(figurefilepath, final_centroids, elec_in_brain):
    """
    Construct PCA plots to visually show deviations in predicted centroids
    and manually labeled centroids.

    Parameters
    ----------
        figurefilepath: str
            File path to save the image to

        final_centroids: dict()
            Dictionary of the centroids for each contact along an electrode
            found by clustering algorithm.

        elec_in_brain: dict()
            Dictionary of manually labeled centroids
    """
    if not os.path.exists(os.path.dirname(figurefilepath)):
        os.mkdir(os.path.dirname(figurefilepath))

    numelectrodes = len(final_centroids.keys())
    print(f"Plotting results for {numelectrodes} electrodes.")

    sns.set(font_scale=1.1)
    fig, axs = plt.subplots(int(np.ceil(numelectrodes / 2)), 2, figsize=(20, 20))
    axs = axs.flatten()

    grouped_chans_labeled = {}
    for label, contact in elec_in_brain.items():
        if label[1] == "'":
            electrode_name = label[:2]  # If electrode name has a '
        else:
            electrode_name = label[0]  # If electrode name has no '
        if electrode_name not in grouped_chans_labeled.keys():
            grouped_chans_labeled[electrode_name] = {}
        grouped_chans_labeled[electrode_name][label] = elec_in_brain[label]
    validation = {}
    for electrode in grouped_chans_labeled:
        validation[electrode] = []
        for channel in grouped_chans_labeled[electrode]:
            validation[electrode].append(grouped_chans_labeled[electrode][channel])

    for i, electrode in enumerate(final_centroids):
        # plot stuff
        pca = PCA()
        pca.fit(np.array(list(final_centroids[electrode].values())))

        centroids_pca = pca.transform(
            np.array(list(final_centroids[electrode].values()))
        )
        centroids_new = pca.inverse_transform(centroids_pca)
        axs[i].scatter(
            centroids_new[:, 0], np.zeros_like(centroids_new[:, 0]), label="observed"
        )
        axs[i].set_title(f"PCA Visualization of Electrode {electrode}")
        pca.fit(validation[electrode])

        expected_centroids_pca = pca.transform(validation[electrode])
        expected_centroids_new = pca.inverse_transform(expected_centroids_pca)
        axs[i].scatter(
            expected_centroids_new[:, 0],
            np.zeros_like(expected_centroids_new[:, 0]),
            marker="x",
            c="r",
            label="expected",
        )
        axs[i].set(
            title=f"PCA Validation of Centroids (Electrode: {electrode})",
            xlabel=f"PC Coordinates in Voxels along Electrode {electrode}",
            ylim=[-0.005, 0.005],
        )
        axs[i].legend()
        for j, chan in enumerate(final_centroids[electrode]):
            axs[i].annotate(chan, (centroids_new[j, 0], 0.0005), size=8.5)

        # set plot paramsx
    fig.tight_layout()

    plt.savefig(figurefilepath, box_inches="tight")

    return fig, axs


def l2_error(figurefilepath, final_centroids, elec_in_brain):
    """
    Function that computes the Euclidean distance (L2) between centroids
    computed by algorithm and the validation data, which are centroids
    manually inputted by user.

    Parameters
    ----------
    final_centroids: dict()
        dictionary of properly labeled centroids grouped by electrode
        in CT voxels.
    elec_in_brain: dict()
        electrode coordinates in CT voxels

    Returns
    -------
        Error in each channel (accurate prediction is between 0-3).
    """
    numelectrodes = len(final_centroids.keys())

    errors_per_channel = {}
    for electrode in final_centroids:
        errors_per_channel[electrode] = {}
        for channel in final_centroids[electrode]:
            observed = final_centroids[electrode][
                channel
            ]  # The coordinates detected by clustering algorithm
            if channel in elec_in_brain:
                expected = elec_in_brain[channel]  # Manually labeled point
                abs_error = npl.norm(observed - expected)
                errors_per_channel[electrode][channel] = abs_error
            else:
                print(
                    f"Channel {channel} from prediction not in brain - saving error as NaN."
                )
                errors_per_channel[electrode][channel] = np.nan

    sns.set(font_scale=1.1)
    fig, axs = plt.subplots(int(np.ceil(numelectrodes / 2.0)), 2, figsize=(15, 15))
    axs = axs.flatten()

    ymin, ymax = 0, 20
    for i, electrode in enumerate(errors_per_channel.keys()):
        y_pos = np.arange(len(errors_per_channel[electrode]))
        axs[i].bar(
            y_pos,
            list(errors_per_channel[electrode].values()),
            align="center",
            alpha=0.9,
        )
        axs[i].set(
            title="Abs. Error By Channel in Electrode {electrode} After Filling Gaps",
            xlabel="Channel",
            xticks=y_pos,
            xticklabels=list(final_centroids[electrode].keys()),
            ylabel="Distance",
            ylim=[ymin, ymax],
        )
        axs[i].set_xlabel("Channel")
        axs[i].set_ylabel("Distance")
        axs[i].set_ylim([ymin, ymax])
    fig.tight_layout()

    plt.savefig(figurefilepath, box_inches="tight")

    return errors_per_channel, fig, axs


def load_data(elecfile):
    """
    Load each brain image scan as a NiBabel image object.

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
    else:
        matreader = MatReader()
        data = matreader.loadmat(elecfile)
        elec_coords_mm = data["data"]

    print("\n\nFinished loading electrode locations data: ")
    print(list(elec_coords_mm.keys()))
    return elec_coords_mm


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "final_centroids_voxels", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "final_centroids_xyz", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "orgfinal_centroids_voxels", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "orgfinal_centroids_xyz", help="The CT image volume in its original space."
    )
    parser.add_argument("ctfile", help="The CT img file")
    parser.add_argument("brainmaskfile", help="The Brain mask file")
    parser.add_argument("fsdir", help="The freesurfer output diretroy.")
    parser.add_argument("summary_pcaplots_filepath")
    parser.add_argument("summary_l2errorplots_filepath")
    parser.add_argument("original_centroids_xyz")
    args = parser.parse_args()

    # extract arguments from parser
    final_centroids_voxels_file = args.final_centroids_voxels
    final_centroids_xyz_file = args.final_centroids_xyz
    orgfinal_centroids_voxels_file = args.orgfinal_centroids_voxels
    orgfinal_centroids_xyz_file = args.orgfinal_centroids_xyz
    ctfile = args.ctfile
    brainmaskfile = args.brainmaskfile
    fsdir = args.fsdir
    pcafig_filepath = args.summary_pcaplots_filepath
    l2fig_filepath = args.summary_l2errorplots_filepath

    original_centroids_xyz = args.original_centroids_xyz

    # Load image data and get the image volumes
    bm_ct_img = nb.load(brainmaskfile)
    ct_img = nb.load(ctfile)
    ct_data = ct_img.get_fdata()
    bm_ct_data = bm_ct_img.get_fdata()

    # load elecs data
    print("LOADING ELECS DATA FROM: ", final_centroids_xyz_file)
    elec_finalcoords_mm = load_data(orgfinal_centroids_xyz_file)

    print("LOADING ORIG ELECS DATA FROM: ", final_centroids_xyz_file)
    elec_coords_mm = load_data(original_centroids_xyz)

    # define pipeline objects to run algorithm
    maskpipe = MaskVolume()

    # Apply masking
    brainmasked_ct_data = maskpipe.apply_mask(bm_ct_data, ct_data)
    brainmasked_ct_img = nb.Nifti1Image(brainmasked_ct_data, bm_ct_img.affine)

    # Filtering out electrodes not within brainmask
    elec_in_brain = maskpipe.filter_electrodes_bm(elec_coords_mm, brainmasked_ct_img)

    """   SUMMARY PLOTS  """
    # create generated figures to check
    fig_pca, axs_pca = summary_PCA_plots(
        pcafig_filepath, elec_finalcoords_mm, elec_in_brain
    )
    l2_errors, fig_l2, axs_l2 = l2_error(
        l2fig_filepath, elec_finalcoords_mm, elec_in_brain
    )

    """   EUCLIDEAN DISTANCE ERROR CHECK """
    l2_errorstxt_path = os.path.join(
        os.path.dirname(pcafig_filepath), "euclidean_distance_errors.txt"
    )
    with open(l2_errorstxt_path, "w") as f:
        for elec in l2_errors:
            for chan in l2_errors[elec]:
                err = l2_errors[elec][chan]
                f.write("%s %.6f\n" % (chan, err))

    # Sum over all individual errors
    total_err = 0
    for elec in l2_errors:
        for chan in l2_errors[elec]:
            if not math.isnan(l2_errors[elec][chan]):
                total_err += l2_errors[elec][chan]

    print(f"Total L2 error in analysis: {total_err}")
