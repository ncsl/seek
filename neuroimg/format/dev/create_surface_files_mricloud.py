import sys

import numpy as np

sys.path.append("../../../")


def read_mricloud_lut(lut_filepath):
    #
    pass


def generate_surfaces():
    pass
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
