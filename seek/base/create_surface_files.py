# -*- coding: utf-8 -*-

import argparse
import os
import os.path
import sys
import time

import numpy as np

sys.path.append("../../")

from seek.base.objects.baseneuroimage import RegionIndexMapping
from seek.base.objects.surface_object import GetSurface
from seek.base.objects.dataset import StructuralDataset

"""
Main creation file
"""


def create_surface_main(
    cort_surf_direc: os.PathLike,
    label_direc: os.PathLike,
    subcort_surf_direc: os.PathLike,
    source_lut: os.PathLike,
    target_lut: os.PathLike,
    struct_zip_file: os.PathLike,
    out_surfaces_dir: os.PathLike = None,
    include_unknown: bool = False,
):
    """
    Create a zipped file with all relevant files for interacting with a Surface of the brain.

    Parameters
    ----------
    cort_surf_direc: Directory that should contain:
                       rh.pial
                       lh.pial
    label_direc: Directory that should contain:
                   rh.aparc.annot
                   lh.aparc.annot
    subcort_surf_direc: Directory that should contain:
                          aseg_<NUM>.srf
                       for each <NUM> in SUBCORTICAL_REG_INDS

    source_lut: File with the color look-up table used for the original parcellation

    target_lut: File with the color look-up table used for the connectome generation

    struct_zip_file: zip file containing TVB structural dataset to be created

    out_surfaces_dir: directory where to put the surfaces and region mappings in TVB format

    include_unknown: include also unknown regions in the connectome
    """
    print("starting...")
    tic = time.time()

    surface_worker = GetSurface

    # create the region index mapping into the target lut space
    region_index_mapping = RegionIndexMapping(source_lut, target_lut)

    # compute the subcortical and cortical surfaces
    surf_subcort = surface_worker.get_subcortical_surfaces(
        subcort_surf_direc, region_index_mapping
    )
    surf_cort = surface_worker.get_cortical_surfaces(
        cort_surf_direc, label_direc, region_index_mapping
    )

    # compute the region parameters for subcort and cort - regions, areas, orientations, centers
    region_params_subcort = surface_worker.compute_region_params(surf_subcort, True)
    region_params_cort = surface_worker.compute_region_params(surf_cort, False)

    nregions = max(region_index_mapping.trg_table.inds) + 1
    orientations = np.zeros((nregions, 3))
    areas = np.zeros(nregions)
    centers = np.zeros((nregions, 3))
    cortical = np.zeros(nregions, dtype=bool)

    # loop through subcortical and cortical regions
    for region_params, is_cortical in [
        (region_params_subcort, False),
        (region_params_cort, True),
    ]:
        regions, reg_areas, reg_orientations, reg_centers = region_params
        orientations[regions, :] = reg_orientations
        areas[regions] = reg_areas
        centers[regions, :] = reg_centers
        cortical[regions] = is_cortical

    if not include_unknown:
        # Remove the region from orientations, areas and centers
        indices = list(range(0, region_index_mapping.unknown_ind)) + list(
            range(region_index_mapping.unknown_ind + 1, nregions)
        )

        names = region_index_mapping.trg_table.names[indices]
        orientations = orientations[indices]
        areas = areas[indices]
        centers = centers[indices]
        cortical = cortical[indices]

        remap_dict = {
            ind: ind if ind < region_index_mapping.unknown_ind else ind - 1
            for ind in range(nregions)
        }
        remap_dict[region_index_mapping.unknown_ind] = -1
        surf_subcort.remap(remap_dict)
        surf_cort.remap(remap_dict)

    else:
        # extract the names of each region index from target lut
        names = region_index_mapping.trg_table.names

    dataset = StructuralDataset(orientations, areas, centers, cortical, names)
    dataset.save_to_txt_zip(struct_zip_file)

    if out_surfaces_dir:
        surf_subcort.save_region_mapping_txt(
            os.path.join(out_surfaces_dir, "region_mapping_subcort.txt")
        )
        surf_subcort.save_surf_zip(
            os.path.join(out_surfaces_dir, "surface_subcort.zip")
        )
        surf_cort.save_region_mapping_txt(
            os.path.join(out_surfaces_dir, "region_mapping_cort.txt")
        )
        surf_cort.save_surf_zip(os.path.join(out_surfaces_dir, "surface_cort.zip"))

    print("complete in %0.2fs", time.time() - tic)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("subject_dir", help="The result directory of the dataset")
    parser.add_argument(
        "subcort_dir", help="The directory with the subcortical surfaces files."
    )
    parser.add_argument(
        "source_lut", help="Figure directory to save resulting figures."
    )
    parser.add_argument("target_lut", help="The output model computed filename.")
    parser.add_argument("struct_zip_file", help="The output meta data filename")
    parser.add_argument(
        "out_surfaces_dir", help="The output path for the figure name to be saved."
    )
    args = parser.parse_args()
    # extract arguments from parser
    subject_dir = args.subject_dir
    source_lut = args.source_lut
    target_lut = args.target_lut
    struct_zip_file = args.struct_zip_file
    out_surfaces_dir = args.out_surfaces_dir
    subcort_surf_direc = args.subcort_dir

    # directories to the cortical and subcortical data
    cort_surf_direc = os.path.join(subject_dir, "surf")
    # subcort_surf_direc = os.path.join(subject_dir, "aseg2srf")
    label_direc = os.path.join(subject_dir, "label")
    include_unknown = True

    create_surface_main(
        cort_surf_direc,
        label_direc,
        subcort_surf_direc,
        source_lut,
        target_lut,
        struct_zip_file,
        out_surfaces_dir,
        include_unknown,
    )
