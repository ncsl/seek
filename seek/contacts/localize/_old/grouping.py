import numpy as np
import numpy.linalg as npl
from skimage import measure


class Cluster:
    """Class of functions for clustering contacts into electrodes."""

    @classmethod
    def cluster_with_threshold(
        self, brainmasked_CT_arr: np.ndarray, threshold: float = 0.630
    ):
        """
        Apply a thresholding based algorithm and then running connected clustering using skimage.

        This will then return clustered_voxels of voxels that belong to distinct contact groups.
        http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label

        Parameters
        ----------
            brainmasked_CT_arr: ndarray
                3-dimensional brain-masked CT scan image.

            threshold: float
                Brightness threshold to use for binarizing the input image.

        Returns
        -------
            clustered_voxels: dict()
                Dictionary that maps cluster ID's to voxels belonging to that
                cluster.

            numobj: int
                Number of clustered_voxels found in image for specified threshold.
        """
        if threshold < 0 or threshold > 1.0:  # pragma: no cover
            raise ValueError(
                "Threshold for clustering should be between 0 and 1, since "
                "threshold is based on the normalized voxel brightness (i.e. "
                f"all voxel values are normalized by 255. You passed in {threshold}."
            )
        if brainmasked_CT_arr.ndim != 3:  # pragma: no cover
            raise RuntimeError(
                "Brain masked CT array should be 3 dimensional? "
                f"The image data you passed in has shape: {brainmasked_CT_arr.shape}"
            )

        # initialize clustered_voxels as a dictionary by ID
        clusters = {}

        print("Threshold \t Num Clusters Found")

        # Label voxels with ID corresponding to the cluster it belongs to
        cluster_labels, numobj = measure.label(
            brainmasked_CT_arr / 255 > threshold,
            background=0,
            return_num=True,
            connectivity=2,
        )

        print(f"{threshold} \t\t {numobj}")

        # Filter out all zero-valued voxels in cluster_labels (pixels)
        nonzeros = np.nonzero(cluster_labels)
        nonzero_voxs = np.array(list(zip(*nonzeros)))

        # Group voxels by their cluster ID
        for vox in nonzero_voxs:
            clusters.setdefault(cluster_labels[tuple(vox)], []).append(vox)

        assert len(clusters) == numobj
        return clusters, numobj


class CylindricalGroup:
    """Class of static functions for cylindrical grouping."""

    @classmethod
    def _is_point_in_cylinder(self, pt1, pt2, radius: int, q):
        """
        Test whether a point q lies within a cylinder.

        With points pt1 and pt2 that define the axis of the cylinder and a
        specified radius r. Used formulas provided here:

        https://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml

        Parameters
        ----------
            pt1: ndarray
                first point to bound the cylinder. (N, 1)

            pt2: ndarray
                second point to bound the cylinder. (N, 1)

            radius: int
                the radius of the cylinder.

            q: ndarray
                the point to test whether it lies within the cylinder.

        Returns
        -------
            True if q lies in the cylinder of specified radius formed by pt1
            and pt2, False otherwise.
        """
        # convert to arrays
        pt1 = np.array(pt1)
        pt2 = np.array(pt2)
        q = np.array(q)

        if any([x.shape != pt1.shape for x in [pt2, q]]):  # pragma: no cover
            print(pt1, pt2, q)
            raise RuntimeError(
                "All points that are checked in cylinder "
                "should have the same shape. "
                f"The shapes passed in were: {pt1.shape}, {pt2.shape}, {q.shape}"
            )

        vec = pt2 - pt1
        length_sq = npl.norm(vec) ** 2
        radius_sq = radius ** 2
        testpt = q - pt1  # set pt1 as origin
        dot = np.dot(testpt, vec)

        # Check if point is within caps of the cylinder
        if dot >= 0 and dot <= length_sq:
            # Distance squared to the cylinder axis of the cylinder:
            dist_sq = np.dot(testpt, testpt) - (dot * dot / length_sq)
            if dist_sq <= radius_sq:
                return True

        return False

    @classmethod
    def apply_cylinder_filter(
        self, sparse_labeled_contacts_vox, clustered_voxels, radius
    ):
        """
        Apply a cylindrical filtering on the raw threshold-based clustering.

        Generates bounding cylinders given a sparse
        collection of contact coordinates.

        Parameters
        ----------
            sparse_labeled_contacts_vox: dict(str: dict(str: ndarray))
                Dictionary of contacts grouped by electrode. An electrode name
                maps to a dictionary of contact labels and corresponding
                coordinates. The dictionary is in sorted order based on these
                labels.

            clustered_voxels: dict()
                Dictionary that maps cluster ID's to voxels belonging to that
                cluster.

            radius: int
                Radius with which to form cylindrical boundaries.

        Returns
        -------
            voxel_clusters: dict(str: dict(str: ndarray))
                Dictionary of clusters sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.
        """
        # each electrode corresponds now to a separate cylinder
        electrode_voxel_clusters = {}

        # go through each electrode apply cylinder filter based on entry/exit point
        for (
            elec,
            (entry_point_ch, exit_point_ch),
        ) in sparse_labeled_contacts_vox.items():
            entry_point_xyz = sparse_labeled_contacts_vox[elec][entry_point_ch]
            exit_point_xyz = sparse_labeled_contacts_vox[elec][exit_point_ch]
            # remove all voxels in clusters that are not within cylinder
            for cluster_id, cluster_voxels in clustered_voxels.items():
                # Obtain all points within the electrode endpoints pt1, pt2
                contained_points = [
                    point
                    for point in cluster_voxels
                    if self._is_point_in_cylinder(
                        entry_point_xyz, exit_point_xyz, radius, point
                    )
                ]

                if not contained_points:
                    continue
                cluster_voxels = np.array(cluster_voxels)
                electrode_voxel_clusters.setdefault(elec, {})[
                    cluster_id
                ] = cluster_voxels

        return electrode_voxel_clusters
