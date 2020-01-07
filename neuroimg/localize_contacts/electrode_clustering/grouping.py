import numpy as np
import numpy.linalg as npl
from skimage import measure


class Cluster:
    @classmethod
    def find_clusters(self, maskedCT, threshold=0.630):
        """
        Function to apply a thresholding based algorithm and then running
        connected clustering using skimage. This will then return clusters
        of voxels that belong to distinct contact groups.
        http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label

        Parameters
        ----------
            maskedCT: ndarray
                3-dimensional brain-masked CT scan image.

            threshold: float
                Brightness threshold to use for binarizing the input image.

        Returns
        -------
            clusters: dict()
                Dictionary that maps cluster ID's to voxels belonging to that
                cluster.

            numobj: int
                Number of clusters found in image for specified threshold.
        """
        clusters = {}

        print("Threshold \t Num Clusters Found")

        # Label voxels with ID corresponding to the cluster it belongs to
        cluster_labels, numobj = measure.label(
            maskedCT / 255 > threshold,
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

        return clusters, numobj


class CylindricalGroup:
    @classmethod
    def _point_in_cylinder(self, pt1, pt2, radius, q):
        """
        Tests whether a point q lies within a cylinder with points
        pt1 and pt2 that define the axis of the cylinder and a
        specified radius. Used formulas provided here:
        https://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml

        Parameters
        ----------
            pt1: ndarray
                first point to bound the cylinder.

            pt2: ndarray
                second point to bound the cylinder.

            radius: int
                the radius of the cylinder.

            q: ndarray
                the point to test whether it lies within the cylinder.

        Returns
        -------
            True if q lies in the cylinder of specified radius formed by pt1
            and pt2, False otherwise.
        """
        pt1 = np.array(pt1)
        pt2 = np.array(pt2)
        q = np.array(q)

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
    def cylinder_filter(self, labeled_contacts, clusters, radius):
        """
        Applies a cylindrical filtering on the raw threshold-based
        clustering by generating bounding cylinders given a sparse
        collection of contact coordinates.

        Parameters
        ----------
            labeled_contacts: dict(str: dict(str: ndarray))
                Dictionary of contacts grouped by electrode. An electrode name
                maps to a dictionary of contact labels and corresponding
                coordinates. The dictionary is in sorted order based on these
                labels.

            clusters: dict()
                Dictionary that maps cluster ID's to voxels belonging to that
                cluster.

            radius: int
                Radius with which to form cylindrical boundaries.

        Returns
        -------
            bound_clusters: dict(str: dict(str: ndarray))
                Dictionary of clusters sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

            sparse_labeled_contacts: dict(str: dict: ndarray)
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.
        """
        # Perform cylindrical filtering on the points detected by
        # threshold-based clustering
        sparse_labeled_contacts = {}
        for elec in labeled_contacts:

            if len(labeled_contacts[elec]) < 2:
                raise ValueError(
                    f"Too few labeled contacts for electrode {elec}."
                    f"Please provide two contacts per electrode."
                )
            chans = list(labeled_contacts[elec].items())
            first, last = chans[0], chans[-1]
            sparse_labeled_contacts[elec] = dict([first, last])

        if not sparse_labeled_contacts:
            return None

        bound_clusters = {}
        for elec in sparse_labeled_contacts:
            pt1, pt2 = list(sparse_labeled_contacts[elec].values())
            for cluster_id in clusters:
                cluster_pts = clusters[cluster_id]
                # Obtain all points within the electrode endpoints pt1, pt2
                contained_points = [
                    point
                    for point in cluster_pts
                    if self._point_in_cylinder(pt1, pt2, radius, point)
                ]

                if not contained_points:
                    continue
                cluster_pts = np.array(cluster_pts)
                bound_clusters.setdefault(elec, {})[cluster_id] = cluster_pts

        return bound_clusters, sparse_labeled_contacts
