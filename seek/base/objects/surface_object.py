import os
import os.path
import tempfile
from typing import List
from zipfile import ZipFile

import numpy as np

from seek.base.objects.baseneuroimage import (
    Hemisphere,
    RegionIndexMapping,
)
from seek.utils import (
    pial_to_verts_and_triangs,
    read_cortical_region_mapping,
)

SUBCORTICAL_REG_INDS = [
    8,
    10,
    11,
    12,
    13,
    16,
    17,
    18,
    26,
    47,
    49,
    50,
    51,
    52,
    53,
    54,
    58,
]


class Surface:
    """
    Class wrapper for the surface data used for being able to reconstruct the surface geometry.

    Takes in vertices, triangles and region_mappings to compute:
    1. vertex triangles
    2. vertex normals
    3. triangle areas
    4. triangle normals
    5. triangle angles
    """

    def __init__(
        self, vertices: np.array, triangles: np.array, region_mapping: np.array
    ):
        assert vertices.ndim == 2
        assert triangles.ndim == 2
        assert region_mapping.ndim == 1

        assert vertices.shape[1] == 3
        assert triangles.shape[1] == 3
        assert region_mapping.shape[0] == vertices.shape[0]

        self.vertices = vertices
        self.triangles = triangles
        self.region_mapping = region_mapping

        self.nverts = self.vertices.shape[0]
        self.ntriangs = self.triangles.shape[0]

        # compute vertex and edges
        self.vertex_triangles = Surface.compute_vertex_triangles(
            self.nverts, self.ntriangs, self.triangles
        )
        self.triangle_normals = Surface.compute_triangle_normals(
            self.triangles, self.vertices
        )
        self.triangle_angles = Surface.compute_triangle_angles(
            self.vertices, self.ntriangs, self.triangles
        )
        self.vertex_normals = Surface.compute_vertex_normals(
            self.nverts,
            self.vertices,
            self.vertex_triangles,
            self.triangles,
            self.triangle_angles,
            self.triangle_normals,
        )
        self.triangle_areas = Surface.compute_triangle_areas(
            self.vertices, self.triangles
        )

    def save_surf_zip(self, filename):
        """
        Create a zipped file with the necessary .txt files to recreate the brain surfaces.

        Parameters
        ----------
        filename :

        """
        tmpdir = tempfile.TemporaryDirectory()

        # save vertices, triangles, and vertex normals
        file_vertices = os.path.join(tmpdir.name, "vertices.txt")
        file_triangles = os.path.join(tmpdir.name, "triangles.txt")
        file_normals = os.path.join(tmpdir.name, "normals.txt")

        np.savetxt(file_vertices, self.vertices, fmt="%.6f %.6f %.6f")
        np.savetxt(file_triangles, self.triangles, fmt="%d %d %d")
        np.savetxt(file_normals, self.vertex_normals, fmt="%.6f %.6f %.6f")

        with ZipFile(filename, "w") as zip_file:
            zip_file.write(file_vertices, os.path.basename(file_vertices))
            zip_file.write(file_triangles, os.path.basename(file_triangles))
            zip_file.write(file_normals, os.path.basename(file_normals))

    def remap(self, remap_dict):
        """
        Remap region mapping to another lookup table.

        Parameters
        ----------
        remap_dict : The other lookup table

        """
        for old_ind, new_ind in remap_dict.items():
            self.region_mapping[self.region_mapping == old_ind] = new_ind

    def save_region_mapping_txt(self, filename):
        """
        Save the region mapping as its own file.

        Parameters
        ----------
        filename : (string) a filename.txt

        """
        np.savetxt(filename, self.region_mapping, fmt="%d")

    @staticmethod
    def compute_vertex_triangles(number_of_vertices, number_of_triangles, triangles):
        """
        Compute the vertex triangles.

        :param number_of_vertices:
        :param number_of_triangles:
        :param triangles:
        :return:
        """
        vertex_triangles = [[] for _ in range(number_of_vertices)]
        for k in range(number_of_triangles):
            vertex_triangles[triangles[k, 0]].append(k)
            vertex_triangles[triangles[k, 1]].append(k)
            vertex_triangles[triangles[k, 2]].append(k)
        return vertex_triangles

    @staticmethod
    def compute_triangle_normals(triangles, vertices):
        """
        Calculate triangle normals.

        :param triangles:
        :param vertices:
        :return:
        """
        tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
        tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
        tri_norm = np.cross(tri_u, tri_v)

        try:
            triangle_normals = (
                tri_norm / np.sqrt(np.sum(tri_norm ** 2, axis=1))[:, np.newaxis]
            )
        except FloatingPointError:
            # TODO: NaN generation would stop execution, however for normals this case could maybe be
            #  handled in a better way.
            triangle_normals = tri_norm
        return triangle_normals

    @staticmethod
    def compute_triangle_angles(vertices, number_of_triangles, triangles):
        """
        Calculate the inner angles of all the triangles which make up a surface.

        TODO: Should be possible with arrays, ie not nested loops...
        A short profile indicates this function takes 95% of the time to compute normals
        (this was a direct translation of some old matlab code)

        :param vertices:
        :param number_of_triangles:
        :param triangles:
        :return:
        """
        verts = vertices
        # TODO: Should be possible with arrays, ie not nested loops...
        # A short profile indicates this function takes 95% of the time to compute normals
        # (this was a direct translation of some old matlab code)
        angles = np.zeros((number_of_triangles, 3))
        for tt in range(number_of_triangles):
            triangle = triangles[tt, :]
            for ta in range(3):
                ang = np.roll(triangle, -ta)
                angles[tt, ta] = np.arccos(
                    np.dot(
                        (verts[ang[1], :] - verts[ang[0], :])
                        / np.sqrt(
                            np.sum((verts[ang[1], :] - verts[ang[0], :]) ** 2, axis=0)
                        ),
                        (verts[ang[2], :] - verts[ang[0], :])
                        / np.sqrt(
                            np.sum((verts[ang[2], :] - verts[ang[0], :]) ** 2, axis=0)
                        ),
                    )
                )

        return angles

    @staticmethod
    def compute_vertex_normals(
        number_of_vertices,
        vertices,
        vertex_triangles,
        triangles,
        triangle_angles,
        triangle_normals,
    ):
        """
        Estimates vertex normals, based on triangle normals weighted by the angle they subtend at each vertex.

        :param number_of_vertices: (int) the number of vertices to analyze
        :param vertices:
        :param vertex_triangles:
        :param triangles:
        :param triangle_angles:
        :param triangle_normals:
        :return: vertex normal vectors
        """
        vert_norms = np.zeros((number_of_vertices, 3))
        bad_normal_count = 0
        for k in range(number_of_vertices):
            try:
                tri_list = list(vertex_triangles[k])
                angle_mask = triangles[tri_list, :] == k
                angles = triangle_angles[tri_list, :]
                angles = angles[angle_mask][:, np.newaxis]
                angle_scaling = angles / np.sum(angles, axis=0)
                vert_norms[k, :] = np.mean(
                    angle_scaling * triangle_normals[tri_list, :], axis=0
                )
                # Scale by angle subtended.
                vert_norms[k, :] = vert_norms[k, :] / np.sqrt(
                    np.sum(vert_norms[k, :] ** 2, axis=0)
                )
                # Normalise to unit vectors.
            except (ValueError, FloatingPointError):
                # If normals are bad, default to position vector
                # A nicer solution would be to detect degenerate triangles and ignore their
                # contribution to the vertex normal
                vert_norms[k, :] = vertices[k] / np.sqrt(vertices[k].dot(vertices[k]))
                bad_normal_count += 1
        if bad_normal_count:
            print(" %d vertices have bad normals" % bad_normal_count)
        return vert_norms

    @staticmethod
    def compute_triangle_areas(vertices, triangles):
        """
        Calculate the area of triangles making up a surface.

        :param vertices:
        :param triangles:
        :return:
        """
        tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
        tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
        tri_norm = np.cross(tri_u, tri_v)
        triangle_areas = np.sqrt(np.sum(tri_norm ** 2, axis=1)) / 2.0
        triangle_areas = triangle_areas[:, np.newaxis]
        return triangle_areas


class GetSurface:
    """Class for surface functions."""

    @staticmethod
    def get_cortical_surfaces(
        cort_surf_direc: os.PathLike,
        label_direc: os.PathLike,
        region_index_mapping: RegionIndexMapping,
    ) -> Surface:
        """
        Compute and format a cortical Surface object.

        :param cort_surf_direc: (string) The directory for the cortical surface with rh.pial and lh.pial
        :param label_direc: (string) The directory with the annotations such as rh.aparc.annot
        :param region_index_mapping: (RegionIndexMapping) The region index mapping
        :return: Surface
        """
        # compute the vertices and triangles for left/right hemisphere regions
        verts_l, triangs_l = pial_to_verts_and_triangs(
            os.path.join(cort_surf_direc, Hemisphere.lh.value + ".pial")
        )
        verts_r, triangs_r = pial_to_verts_and_triangs(
            os.path.join(cort_surf_direc, Hemisphere.rh.value + ".pial")
        )

        # get the cortical region mapping for left/right hemispheres
        region_mapping_l = read_cortical_region_mapping(
            label_direc, Hemisphere.lh, region_index_mapping
        )
        region_mapping_r = read_cortical_region_mapping(
            label_direc, Hemisphere.rh, region_index_mapping
        )

        # construct the surface that is concatenation between the two hemispheres
        surface = GetSurface.merge_surfaces(
            [
                Surface(verts_l, triangs_l, region_mapping_l),
                Surface(verts_r, triangs_r, region_mapping_r),
            ]
        )

        return surface

    @staticmethod
    def get_subcortical_surfaces(
        subcort_surf_direc: os.PathLike, region_index_mapping: RegionIndexMapping
    ) -> Surface:
        """
        Compute and format the subcortical Surface.

        :param subcort_surf_direc: The directory for the subcortical surface.
        :param region_index_mapping: (RegionIndexMapping) The region index mapping for subcortical structs
        :return: Surface
        """
        # store surfaces in a list
        surfaces = []

        # loop through all SUBCORTICAL REGION INDICES
        for fs_idx in SUBCORTICAL_REG_INDS:
            conn_idx = region_index_mapping.source_to_target(fs_idx)
            filename = os.path.join(subcort_surf_direc, "aseg_%03d.srf" % fs_idx)
            with open(filename, "r") as f:
                f.readline()
                nverts, ntriangs = [int(n) for n in f.readline().strip().split(" ")]

            vertices = np.genfromtxt(
                filename,
                dtype=float,
                skip_header=2,
                skip_footer=ntriangs,
                usecols=(0, 1, 2),
            )
            triangles = np.genfromtxt(
                filename, dtype=int, skip_header=2 + nverts, usecols=(0, 1, 2)
            )
            region_mapping = conn_idx * np.ones(nverts, dtype=int)
            surfaces.append(Surface(vertices, triangles, region_mapping))

        surface = GetSurface.merge_surfaces(surfaces)
        return surface

    @staticmethod
    def merge_surfaces(surfaces: List[Surface]) -> Surface:
        """
        Merge a list of Surfaces into one Surface object.

        :param surfaces: (list[Surface]) a list of surfaces
        :return: Surface
        """
        offsets = np.cumsum(
            [0] + [vs.shape[0] for vs in [surf.vertices for surf in surfaces]][:-1]
        )
        vertices = np.vstack([surf.vertices for surf in surfaces])
        triangles = np.vstack(
            [
                ts + offset
                for ts, offset in zip([surf.triangles for surf in surfaces], offsets)
            ]
        )
        region_mappings = np.hstack([surf.region_mapping for surf in surfaces])
        return Surface(vertices, triangles, region_mappings)

    @staticmethod
    def compute_region_params(
        surface: Surface, subcortical: bool = False
    ) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
        """
        Compute the parameters required for each region (regions, areas, orientations, centers).

        Parameters
        ----------
        surface : (Surface) the surface object.
        subcortical : (bool) whether or not this surface is subcortical

        Returns
        -------
        regions, areas, orientations, centers : (tuple of np.ndarray) the regions, areas, orientations and the volume centers
        """
        verts, triangs, region_mapping = (
            surface.vertices,
            surface.triangles,
            surface.region_mapping,
        )

        regions = np.unique(region_mapping)
        areas = GetSurface.compute_region_areas(
            regions, surface.triangle_areas, surface.vertex_triangles, region_mapping
        )
        orientations = GetSurface.compute_region_orientations(
            regions, surface.vertex_normals, region_mapping
        )
        centers = GetSurface.compute_region_centers(regions, verts, region_mapping)

        return regions, areas, orientations, centers

    @staticmethod
    def compute_region_orientations(regions, vertex_normals, region_mapping):
        """
        Compute the orientation of given regions from vertex_normals and region mapping.

        Parameters
        ----------
        regions :
        vertex_normals :
        region_mapping :

        Returns
        -------
        average_orientation :

        """
        average_orientation = np.zeros((regions.size, 3))
        # Average orientation of the region
        for i, k in enumerate(regions):
            orient = vertex_normals[region_mapping == k, :]
            if orient.shape[0] > 0:
                avg_orient = np.mean(orient, axis=0)
                average_orientation[i, :] = avg_orient / np.sqrt(
                    np.sum(avg_orient ** 2)
                )

        return average_orientation

    @staticmethod
    def compute_region_areas(regions, triangle_areas, vertex_triangles, region_mapping):
        """
        Compute the areas of given regions.

        Parameters
        ----------
        regions :
        triangle_areas :
        vertex_triangles :
        region_mapping :

        Returns
        -------
        region_surface_area :

        """
        region_surface_area = np.zeros(regions.size)
        avt = np.array(vertex_triangles)
        # NOTE: Slightly overestimates as it counts overlapping border triangles,
        #       but, not really a problem provided triangle-size << region-size.
        for i, k in enumerate(regions):
            regs = map(set, avt[region_mapping == k])
            region_triangles = set.union(*regs)
            if region_triangles:
                region_surface_area[i] = triangle_areas[list(region_triangles)].sum()

        return region_surface_area

    @staticmethod
    def compute_region_centers(regions, vertices, region_mapping):
        """
        Compute region centers from their vertices and region mapping.

        Parameters
        ----------
        regions :
        vertices :
        region_mapping :

        Returns
        -------
        region_centers :
        """
        region_centers = np.zeros((regions.size, 3))
        for i, k in enumerate(regions):
            vert = vertices[region_mapping == k, :]
            if vert.shape[0] > 0:
                region_centers[i, :] = np.mean(vert, axis=0)

        return region_centers
