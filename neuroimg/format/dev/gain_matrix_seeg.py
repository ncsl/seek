#!/usr/bin/env python3

import argparse
import os
import zipfile

import numpy as np

SIGMA = 1.0


def gain_matrix_dipole(vertices: np.ndarray, orientations: np.ndarray, areas: np.ndarray, region_mapping: np.ndarray,
                       nregions: int, sensors: np.ndarray):
    """

    Parameters
    ----------
    vertices             np.ndarray of floats of size n x 3, where n is the number of vertices
    orientations         np.ndarray of floats of size n x 3
    region_mapping       np.ndarray of ints of size n
    sensors              np.ndarray of floats of size m x 3, where m is the number of sensors

    Returns
    -------
    np.ndarray of size m x n

    """

    nverts = vertices.shape[0]
    nsens = sensors.shape[0]

    reg_map_mtx = np.zeros((nverts, nregions), dtype=int)
    for i, region in enumerate(region_mapping):
        if region >= 0:
            reg_map_mtx[i, region] = 1
    # reg_map_mtx[np.arange(region_mapping.size), region_mapping] = 1.0

    gain_mtx_vert = np.zeros((nsens, nverts))
    for sens_ind in range(nsens):
        a = sensors[sens_ind, :] - vertices
        na = np.sqrt(np.sum(a ** 2, axis=1))
        gain_mtx_vert[sens_ind, :] = areas * (np.sum(orientations * a, axis=1) / na ** 3) / (4.0 * np.pi * SIGMA)

    return gain_mtx_vert @ reg_map_mtx


def gain_matrix_inv_square(vertices: np.ndarray, areas: np.ndarray, region_mapping: np.ndarray,
                           nregions: int, sensors: np.ndarray):
    nverts = vertices.shape[0]
    nsens = sensors.shape[0]

    reg_map_mtx = np.zeros((nverts, nregions), dtype=int)
    for i, region in enumerate(region_mapping):
        if region >= 0:
            reg_map_mtx[i, region] = 1

    gain_mtx_vert = np.zeros((nsens, nverts))
    for sens_ind in range(nsens):
        a = sensors[sens_ind, :] - vertices
        na = np.sqrt(np.sum(a ** 2, axis=1))
        gain_mtx_vert[sens_ind, :] = areas / na ** 2

    return gain_mtx_vert @ reg_map_mtx


def compute_triangle_areas(vertices, triangles):
    """Calculates the area of triangles making up a surface."""
    tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
    tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
    tri_norm = np.cross(tri_u, tri_v)
    triangle_areas = np.sqrt(np.sum(tri_norm ** 2, axis=1)) / 2.0
    triangle_areas = triangle_areas[:, np.newaxis]
    return triangle_areas


def compute_vertex_areas(vertices, triangles):
    triangle_areas = compute_triangle_areas(vertices, triangles)
    vertex_areas = np.zeros((vertices.shape[0]))
    for triang, vertices in enumerate(triangles):
        for i in range(3):
            vertex_areas[vertices[i]] += 1. / 3. * triangle_areas[triang]

    return vertex_areas


def read_surf(directory: os.PathLike, use_subcort):
    reg_map_cort = np.genfromtxt((os.path.join(directory, "region_mapping_cort.txt")), dtype=int)
    reg_map_subc = np.genfromtxt((os.path.join(directory, "region_mapping_subcort.txt")), dtype=int)

    with zipfile.ZipFile(os.path.join(directory, "surface_cort.zip")) as zip:
        with zip.open('vertices.txt') as fhandle:
            verts_cort = np.genfromtxt(fhandle)
        with zip.open('normals.txt') as fhandle:
            normals_cort = np.genfromtxt(fhandle)
        with zip.open('triangles.txt') as fhandle:
            triangles_cort = np.genfromtxt(fhandle, dtype=int)

    with zipfile.ZipFile(os.path.join(directory, "surface_subcort.zip")) as zip:
        with zip.open('vertices.txt') as fhandle:
            verts_subc = np.genfromtxt(fhandle)
        with zip.open('normals.txt') as fhandle:
            normals_subc = np.genfromtxt(fhandle)
        with zip.open('triangles.txt') as fhandle:
            triangles_subc = np.genfromtxt(fhandle, dtype=int)

    vert_areas_cort = compute_vertex_areas(verts_cort, triangles_cort)
    vert_areas_subc = compute_vertex_areas(verts_subc, triangles_subc)

    if not use_subcort:
        return (verts_cort, normals_cort, vert_areas_cort, reg_map_cort)
    else:
        verts = np.concatenate((verts_cort, verts_subc))
        normals = np.concatenate((normals_cort, normals_subc))
        areas = np.concatenate((vert_areas_cort, vert_areas_subc))
        regmap = np.concatenate((reg_map_cort, reg_map_subc))

        return (verts, normals, areas, regmap)


def read_regions(zip_name: os.PathLike, use_subcort):
    with zipfile.ZipFile(zip_name) as zip:
        with zip.open('centres.txt') as fhandle:
            verts = np.genfromtxt(fhandle, usecols=[1, 2, 3])
        with zip.open('areas.txt') as fhandle:
            areas = np.genfromtxt(fhandle)
        with zip.open('average_orientations.txt') as fhandle:
            normals = np.genfromtxt(fhandle)
        with zip.open('cortical.txt') as fhandle:
            cortical = np.genfromtxt(fhandle, dtype=int).astype(bool)
    regmap = np.arange(0, verts.shape[0])

    if not use_subcort:
        return (verts[cortical], normals[cortical], areas[cortical], regmap[cortical])
    else:
        return (verts, normals, areas, regmap)


def get_nregions(zip_name):
    with zipfile.ZipFile(zip_name) as zip:
        with zip.open('centres.txt') as fhandle:
            num_lines = sum(1 for line in fhandle.readlines() if line.strip())
    return num_lines


def main():
    parser = argparse.ArgumentParser(description="Generate SEEG gain matrix.")

    # Defaults are not given on purpose to force the user to think about what is needed.
    parser.add_argument('--mode', type=str, choices=['surface', 'region'], required=True)
    parser.add_argument('--formula', type=str, choices=['dipole', 'inv_square'], required=True)
    parser.add_argument('--surf_dir', help="Directory with surfaces and region mapping. Required if mode is 'surface'.")

    use_subcort_parser = parser.add_mutually_exclusive_group(required=True)
    use_subcort_parser.add_argument('--use_subcort', dest='use_subcort', action='store_true')
    use_subcort_parser.add_argument('--no_use_subcort', dest='use_subcort', action='store_false')

    parser.add_argument('tvb_zipfile', help="Path to the TVB zipfile.")
    parser.add_argument('sensors_file', help="Path to the sensors file.")
    parser.add_argument('gain_matrix', help="Path to the gain matrix in numpy format to be generated.")

    args = parser.parse_args()
    if args.mode == 'surface' and args.surf_dir is None:
        parser.error("--surf_dir is required if mode is 'surface'")

    nregions = get_nregions(args.tvb_zipfile)

    if args.mode == 'surface':
        verts, normals, areas, regmap = read_surf(args.surf_dir, args.use_subcort)
    elif args.mode == 'region':
        verts, normals, areas, regmap = read_regions(args.tvb_zipfile, args.use_subcort)

    sensors_pos = np.genfromtxt(args.sensors_file, usecols=[1, 2, 3])

    if args.formula == 'dipole':
        gain_mtx = gain_matrix_dipole(verts, normals, areas, regmap, nregions, sensors_pos)
    elif args.formula == 'inv_square':
        gain_mtx = gain_matrix_inv_square(verts, areas, regmap, nregions, sensors_pos)

    np.savetxt(args.gain_matrix, gain_mtx)


if __name__ == '__main__':
    main()
