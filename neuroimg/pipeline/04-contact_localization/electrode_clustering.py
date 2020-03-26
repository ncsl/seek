import argparse
import sys
from pathlib import Path

import nibabel as nb
import numpy as np

from .utils.io import save_organized_elecdict_asmat, save_organized_elecdict_astsv, load_elecs_data

sys.path.append("../../../")

from neuroimg.localize_contacts.electrode_clustering.mask import MaskVolume
from neuroimg.localize_contacts.electrode_clustering.grouping import (
    Cluster,
    CylindricalGroup,
)
from neuroimg.localize_contacts.electrode_clustering.postprocess import PostProcessor


def load_data(ct_scan, brainmask_ct=None):
    """
    Load brain image scans as NiBabel image objects.

    Parameters
    –––-------
        ct_scan: str
            Path to Nifti image file of CT scan.

        brainmask_ct: str
            Path to Nifti image file of corresponding brain mask in CT voxels.

    Returns
    -------
        ct_img: NiBabel image object
            NiBabel image object of CT scan input.

        bm_ct_img: NiBabel image object
            NiBabel image object of brain mask in CT.
    """

    ct_img, bm_ct_img = None, None

    if ct_scan:
        ct_img = nb.load(ct_scan)

    if brainmask_ct:
        bm_ct_img = nb.load(brainmask_ct)

    return ct_img, bm_ct_img


def compute_electrode_voxel_clouds(
        elec_coords_mm, ct_img, bm_ct_img=None, output_bm_fpath=None
):
    """Compute electrode voxel clouds via Grouping algo and apply brainmask."""

    # Define pipeline objects to run algorithm
    maskpipe = MaskVolume()

    # Apply masking
    if bm_ct_img:
        brainmasked_ct_img = maskpipe.apply_mask(ct_img, bm_ct_img)
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        # Filtering out electrodes not within brainmask
        elecvoxels_in_brain = maskpipe.mask_electrodes(
            elec_coords_mm, brainmasked_ct_img
        )

        if output_bm_fpath:
            print(f"Saving binarized CT image volume to: {output_bm_fpath}.")
            brainmasked_ct_img.to_filename(output_bm_fpath)
    else:
        if output_bm_fpath:
            raise RuntimeError(
                "Note that output brainmask filepath was passed to intend "
                "saving the brainmasked CT image, but no brainmask Image was passed in."
            )
        brainmasked_ct_data = ct_img.get_fdata()
        elecvoxels_in_brain = maskpipe.mask_electrodes(elec_coords_mm, ct_img)

    # Get all voxel clouds per electrode in sorted order
    grouped_contacts_vox = maskpipe.group_contacts(elecvoxels_in_brain)
    return brainmasked_ct_data, grouped_contacts_vox


def apply_threshold_clustering(brainmasked_ct_data, grouped_contacts_vox, **kwargs):
    # Define pipeline objects to run algorithm
    clusterpipe = Cluster()
    grouppipe = CylindricalGroup()

    threshold = kwargs.get(
        "threshold", 0.630
    )  # radius (in CT voxels) of cylindrical boundary
    radius = kwargs.get(
        "radius", 4
    )  # Between 0 and 1. Zeroes voxels with value < threshold

    # Runs threshold-based clustering algorithm over brainmasked CT
    # for thresholds between 0.63 and 0.64 with a step of 0.005
    clusters, numobj = np.array(
        clusterpipe.find_clusters(brainmasked_ct_data, threshold=threshold)
    )

    # Cluster by cylindrical boundaries
    bound_clusters, sparse_labeled_contacts = grouppipe.cylinder_filter(
        grouped_contacts_vox, clusters, radius
    )

    return bound_clusters, sparse_labeled_contacts


def apply_postprocessing_on_clusters(
        bound_clusters, sparse_labeled_contacts, ct_affine, **kwargs
):
    gap_tolerance = kwargs.get(
        "gap_tolerance", 13
    )  # maximum distance between two adjacent nodes

    # Define pipeline objects to run algorithm
    postprocesspipe = PostProcessor()

    # Begin postprocessing steps
    processed_clusters = postprocesspipe.process_abnormal_clusters(
        bound_clusters, sparse_labeled_contacts
    )

    # Compute centroids for filling gaps
    final_clusters = postprocesspipe.fill_gaps(
        processed_clusters, gap_tolerance, sparse_labeled_contacts
    )

    # Assign channel labels to the clusters found
    final_labeled_clusters = postprocesspipe.assign_labels(
        final_clusters, sparse_labeled_contacts
    )

    # Correct egregious errors with linearly interpolated points
    final_labeled_clusters = postprocesspipe.bruteforce_correction(
        final_labeled_clusters, sparse_labeled_contacts
    )

    # Compute centroids for each cluster
    final_centroids_vox = postprocesspipe.cluster_2_centroids(final_labeled_clusters)

    # Convert final voxels to xyz coordinates
    final_centroids_xyz = postprocesspipe.vox_2_xyz(final_centroids_vox, ct_affine)

    return final_centroids_vox, final_centroids_xyz


def main(ctimgfile, brainmaskfile, elecinitfile, brainmasked_ct_fpath=None):
    # Hyperparameters for electrode clustering algorithm
    radius = 4  # radius (in CT voxels) of cylindrical boundary
    threshold = 0.630  # Between 0 and 1. Zeroes voxels with value < threshold
    gap_tolerance = 13  # maximum distance between two adjacent nodes

    threshold_kwargs = {
        "radius": radius,
        "threshold": threshold,
    }
    postprocess_kwargs = {"gap_tolerance": gap_tolerance}

    # Load data
    print(f"LOADING ELECS DATA FROM: {elecinitfile}")
    elec_coords_mm = load_elecs_data(elecinitfile)
    ct_img, bm_ct_img = load_data(ctimgfile, brainmaskfile)
    ct_affine = ct_img.affine

    # compute electrode voxel clouds and brain masked CT data
    brainmasked_ct_data, grouped_contacts_vox = compute_electrode_voxel_clouds(
        elec_coords_mm, ct_img, bm_ct_img, output_bm_fpath=brainmasked_ct_fpath
    )

    # apply threshold clustering algorithm to sparsely label contacts within clusters
    bound_clusters, sparse_labeled_contacts = apply_threshold_clustering(
        brainmasked_ct_data, grouped_contacts_vox, **threshold_kwargs
    )

    final_centroids_vox, final_centroids_xyz = apply_postprocessing_on_clusters(
        bound_clusters, sparse_labeled_contacts, ct_affine, **postprocess_kwargs
    )

    return (
        final_centroids_vox,
        final_centroids_xyz,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ct_nifti_img", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "brainmask_native_file", help="Brain mask mapped to the CT image space."
    )
    parser.add_argument(
        "electrode_initialization_file",
        help="The electrode file with contacts localized to 2 points.",
    )
    parser.add_argument(
        "clustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument("clustered_voxels_file", help="the voxels output datafile")
    parser.add_argument(
        "orgclustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument("orgclustered_voxels_file", help="the voxels output datafile")

    parser.add_argument("binarized_ct_volume", help="The binarized CT volume.")
    parser.add_argument("fs_patient_dir", help="The freesurfer output diretroy.")
    args = parser.parse_args()

    # Extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    brainmask_native_file = args.brainmask_native_file
    electrode_initialization_file = args.electrode_initialization_file
    clustered_points_file = args.clustered_points_file
    clustered_voxels_file = args.clustered_voxels_file
    binarized_ct_file = args.binarized_ct_volume
    fs_patient_dir = args.fs_patient_dir

    # Create electrodes directory if not exist
    elecs_dir = Path(Path(fs_patient_dir) / "elecs")
    elecs_dir.mkdir(exist_ok=True)

    # Compute the final centroid voxels, centroid xyzs and the binarized CT image volume
    final_centroids_voxels, final_centroids_xyz = main(
        ct_nifti_img,
        brainmask_native_file,
        electrode_initialization_file,
        binarized_ct_file,
    )

    # Save output files
    print(f"Saving clustered xyz coords to: {clustered_points_file}.")
    print(f"Saving clustered voxels to: {clustered_voxels_file}.")
    # save data into bids sidecar-tsv files
    save_organized_elecdict_astsv(final_centroids_xyz, clustered_points_file)
    # save_organized_elecdict_astsv(final_centroids_voxels, clustered_voxels_file)

    # Save centroids as .mat file with attributes eleclabels
    save_organized_elecdict_asmat(final_centroids_xyz, clustered_points_file)
    save_organized_elecdict_asmat(final_centroids_voxels, clustered_voxels_file)
