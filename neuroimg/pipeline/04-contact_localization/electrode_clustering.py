import argparse
import collections
import re
import sys
import warnings
from pathlib import Path
from pprint import pprint

import nibabel as nb

sys.path.append("../../../")

from neuroimg.base.utils.utils import _contact_numbers_on_electrode
from neuroimg.base.utils.io import (
    save_organized_elecdict_asmat,
    save_organized_elecdict_astsv,
    load_elecs_data,
)

from neuroimg.localize_contacts.electrode_clustering.neuroimage import (
    ClusteredBrainImage,
)
from neuroimg.localize_contacts.electrode_clustering.electrode import Electrodes


def get_entry_exit_contacts(electrodes: Electrodes):
    entry_exit_elec = collections.defaultdict(list)
    for elec in electrodes:
        if len(elec) < 6:
            warnings.warn(
                f"Channels on electrode {elec} contain less than 6 contacts - {len(elec)}. "
                "Were endpoints correctly labeled?"
            )

        # get entry/exit contacts
        entry_ch = elec.get_entry_ch()
        exit_ch = elec.get_exit_ch()
        entry_exit_elec[elec.name] = [entry_ch, exit_ch]

        # remove all but the entry/exit
        for ch in elec.contacts:
            if ch.name not in [entry_ch.name, exit_ch.name]:
                elec.remove_contact(ch.name)

    return entry_exit_elec


def mainv2(
    ctimgfile, brainmaskfile, elecinitfile, brainmasked_ct_fpath=None,
):
    # Hyperparameters for electrode clustering algorithm
    radius = 4  # radius (in CT voxels) of cylindrical boundary
    threshold = 0.630  # Between 0 and 1. Zeroes voxels with value < threshold
    contact_spacing_mm = 3.5  # distance between two adjacent contacts

    # load in nibabel CT and brainmask images
    ct_img = nb.load(ctimgfile)
    brainmask_img = nb.load(brainmaskfile)

    # initialize clustered brain image
    brain = ClusteredBrainImage(ct_img, brainmask_img)
    brain.save_masked_img(brainmasked_ct_fpath)  # save brain-masked CT file

    # load in the channel coordinates in xyz as dictionary
    ch_coords_mm = load_elecs_data(elecinitfile)

    # convert into electrodes
    ch_names = list(ch_coords_mm.keys())
    ch_coords = list(ch_coords_mm.values())
    electrodes = Electrodes(ch_names, ch_coords, coord_type="mm")

    # get the entry/exit electrodes
    entry_exit_elec = get_entry_exit_contacts(electrodes)
    # determine the contact numbering per electrode
    elec_contact_nums = {}
    for elec, (entry_ch, exit_ch) in entry_exit_elec.items():
        elec_contact_nums[elec] = _contact_numbers_on_electrode(
            entry_ch.name, exit_ch.name
        )

    # get sparse electrodes in voxel space
    ch_names = []
    ch_coords = []
    # transform coordinates -> voxel space
    for elec_name, contacts in entry_exit_elec.items():
        for contact in contacts:
            contact.transform_coords(brain.get_masked_img(), coord_type="vox")
            ch_names.append(contact.name)
            ch_coords.append(contact.coord)
            assert contact.coord_type == "vox"
    entry_exit_elec = Electrodes(ch_names, ch_coords, coord_type="vox")

    print("Applying SEEK algo... for electrodes: ", entry_exit_elec)
    print("Contact numbering for each electrode: ", elec_contact_nums.keys())
    # compute voxel clusters in the brain image using a threshold
    voxel_clusters, num_clusters = brain.compute_clusters_with_threshold(
        threshold=threshold
    )

    # feed in entry/exit voxel points per electrode and apply a cylinder filter
    _cylindrical_clusters = {}
    for electrode in entry_exit_elec:
        elec_name = electrode.name
        entry_point_vox = electrode.get_entry_ch().coord
        exit_point_vox = electrode.get_exit_ch().coord
        voxel_clusters_in_cylinder = brain.compute_cylindrical_clusters(
            voxel_clusters, entry_point_vox, exit_point_vox, radius=radius
        )
        _cylindrical_clusters[elec_name] = voxel_clusters_in_cylinder
    voxel_clusters = _cylindrical_clusters
    print("Cylindrical bounded electrode clustered_voxels: ", voxel_clusters.keys())

    # preliminarily label electrode voxel clusters
    labeled_voxel_clusters = {}
    electrodes_with_problem = collections.defaultdict(list)
    for electrode in entry_exit_elec:
        elec_name = electrode.name
        entry_ch = electrode.get_entry_ch()
        exit_ch = electrode.get_exit_ch()
        this_elec_voxels = voxel_clusters[elec_name]

        _this_elec_voxels = brain.assign_sequential_labels(
            this_elec_voxels, entry_ch.name, entry_ch.coord
        )

        # check each labels
        for ch_name in _this_elec_voxels.keys():
            _, ch_num = re.match("^([A-Za-z]+[']?)([0-9]+)$", ch_name).groups()
            if int(ch_num) not in elec_contact_nums[elec]:
                electrodes_with_problem[elec_name].append(ch_num)

        # there were incorrectly grouped contact clusters
        if elec_name in electrodes_with_problem:
            print(
                f"Electrode {elec_name} has incorrectly grouped clusters... "
                f"Attempting to unfuse them."
            )
            print(this_elec_voxels.keys())

            # check for merged clusters at entry and exit for this electrode
            merged_cluster_ids = brain._identify_merged_voxel_clusters(this_elec_voxels)

            print(merged_cluster_ids)
            # if there are merged cluster ids, unfuse them
            this_elec_voxels = brain._unfuse_clusters_on_entry_and_exit(
                this_elec_voxels, merged_cluster_ids, elec_contact_nums[elec]
            )

        # check for oversized clusters and pare them down
        oversized_clusters_ids = brain._identify_skull_voxel_clusters(this_elec_voxels)

        print("Found oversized clusters: ", oversized_clusters_ids)

        # pare them down and resize
        this_elec_voxels = brain._pare_clusters_on_electrode(
            this_elec_voxels, oversized_clusters_ids, qtile=0.5
        )

        # fill in gaps between centroids
        this_elec_voxels = brain.fill_clusters_with_spacing(
            this_elec_voxels,
            entry_ch,
            elec_contact_nums[elec_name],
            contact_spacing_mm=contact_spacing_mm,
        )

        if entry_ch.coord_type == 'vox':
            entry_ch.transform_coords(brain.get_masked_img(), coord_type='mm')
        if exit_ch.coord_type == 'vox':
            exit_ch.transform_coords(brain.get_masked_img(), coord_type='mm')
        this_elec_xyz = collections.defaultdict(list)
        for _cluster_id, voxels in this_elec_voxels.items():
            for coord in voxels:
                this_elec_xyz[_cluster_id].append(brain.map_coordinates(coord, coord_type='mm'))
        # # compute the average / std contact-to-contact spacing
        # import numpy as np
        # dists = []
        # for cluster_id, voxels in this_elec_xyz.items():
        #     curr_centroid = np.mean(voxels, axis=0)
        #     if dists == []:
        #         prev_centroid = np.mean(voxels, axis=0)
        #         dists.append(0)
        #         continue
        #     dists.append(np.linalg.norm(curr_centroid - prev_centroid))
        #     prev_centroid = curr_centroid
        # print("Distribution of contact to contact spacing: ", np.mean(dists), np.std(dists))
        # this_elec_xyz = brain.correct_labeled_clusters(this_elec_xyz,
        #                                                   entry_ch,
        #                                                   exit_ch,
        #                                                   contact_spacing_mm=contact_spacing_mm)
        #

        # apply brute force correction
        # this_elec_xyz = brain.bruteforce_correctionv2(
        #     this_elec_xyz,
        #     entry_ch,
        #     exit_ch,
        #     contact_spacing_mm=contact_spacing_mm, num_contacts=len(elec_contact_nums[electrode.name])
        # )
        # _this_elec_voxels = collections.defaultdict(list)
        # for _cluster_id, coords in this_elec_xyz.items():
        #     for coord in coords:
        #         _this_elec_voxels[_cluster_id].append(brain.map_coordinates(coord, coord_type='vox'))
        # this_elec_voxels = _this_elec_voxels

        if entry_ch.coord_type == 'mm':
            entry_ch.transform_coords(brain.get_masked_img(), coord_type='vox')
        if exit_ch.coord_type == 'mm':
            exit_ch.transform_coords(brain.get_masked_img(), coord_type='vox')
        # assign sequential labels
        this_elec_voxels = brain.assign_sequential_labels(
            this_elec_voxels, entry_ch.name, entry_ch.coord,
        )

        # reset electrode clusters to that specific electrode
        labeled_voxel_clusters[elec_name] = this_elec_voxels

    # Compute centroids for each cluster
    labeled_voxel_centroids = brain.cluster_2_centroids(labeled_voxel_clusters)

    # Convert final voxels to xyz coordinates
    labeled_xyz_centroids = brain.vox_2_xyz(
        labeled_voxel_centroids, brain.get_masked_img().affine
    )

    # keep the end points, since they're already labeled
    if entry_exit_elec.coord_type == 'vox':
        entry_exit_elec.transform_coords(brain.get_masked_img(), coord_type='mm')

    for electrode in entry_exit_elec:
        for contact in electrode.contacts:
            labeled_xyz_centroids[electrode.name][contact.name] = contact.coord

    return labeled_voxel_centroids, labeled_xyz_centroids


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
    final_centroids_voxels, final_centroids_xyz = mainv2(
        ct_nifti_img,
        brainmask_native_file,
        electrode_initialization_file,
        binarized_ct_file,
    )

    # Save output files
    print(f"Saving clustered xyz coords to: {clustered_points_file}.")
    print(f"Saving clustered voxels to: {clustered_voxels_file}.")
    pprint(final_centroids_xyz)

    # save data into bids sidecar-tsv files
    save_organized_elecdict_astsv(final_centroids_xyz, clustered_points_file)
    save_organized_elecdict_astsv(final_centroids_voxels, clustered_voxels_file)

    # Save centroids as .mat file with attributes eleclabels
    save_organized_elecdict_asmat(
        final_centroids_xyz, clustered_points_file.replace(".tsv", ".mat")
    )
    save_organized_elecdict_asmat(
        final_centroids_voxels, clustered_voxels_file.replace(".tsv", ".mat")
    )
