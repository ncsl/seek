import collections
import re
from typing import Dict, List, Tuple

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine
from scipy.optimize import linprog
from scipy.stats import norm
from skimage import measure
from sklearn.cluster import KMeans

from seek.contacts.localize.electrode import Contact


def _multidimdist(point1, point2):
    """Compute distance between two N-dimensional points."""
    # find the difference between the two points, its really the same as below
    deltaVals = [
        point2[dimension] - point1[dimension] for dimension in range(len(point1))
    ]
    runningSquared = 0
    # because the pythagorean theorem works for any dimension we can just use that
    for coord in deltaVals:
        runningSquared += coord ** 2
    return runningSquared ** (1 / 2)


def _compute_vector(point1, point2, unitSphere=False):
    """Compute (unit)-vector between two N-dimensional points."""
    finalVector = [0 for coord in point1]
    for dimension, coord in enumerate(point1):
        # finding total difference for that coordinate(x,y,z...)
        deltaCoOrd = point2[dimension] - coord
        # adding total difference
        finalVector[dimension] = deltaCoOrd

    # setting unitSphere to True will make the vector
    # scaled down to a sphere with a radius one,
    # instead of it's original length
    if unitSphere:
        totalDist = _multidimdist(point1, point2)
        unitVector = []
        for dimen in finalVector:
            unitVector.append(dimen / totalDist)
        return unitVector
    else:
        return finalVector


def _in_hull(points, x):
    """Check if x is inside convex hull of points."""
    points = np.array(points)
    x = np.array(x)

    if len(points) == 1:
        return np.array_equal(x, points)

    n_points = len(points)
    c = np.zeros(n_points)
    A = np.r_[points.T, np.ones((1, n_points))]
    b = np.r_[x, np.ones(1)]
    lp = linprog(c, A_eq=A, b_eq=b)
    return lp.success


class BrainImage:
    """
    Apply a binarized brain mask to an input CT image file.

    Parameters
    ----------
        img: NiBabel Image Object
            CT brain image - NIFTI formatting standard.

        mask_img: NiBabel Image Object
            Corresponding brain mask - NIFTI formatting standard.

    """

    def __init__(self, img: nb.Nifti2Image, mask_img: nb.Nifti2Image):
        self.img = img
        self.mask_img = mask_img

        if img.ndim != 3:  # pragma: no cover
            raise RuntimeError(
                "Image data should be 3 dimensional? "
                f"The image data you passed in has shape: {img.shape}"
            )

        if img.shape != mask_img.shape:  # pragma: no cover
            raise ValueError(
                "Brain image and Mask image shapes need to be the "
                f"same. You passed in {img.shape} and  {mask_img.shape}."
            )

    @property
    def mask_arr(self):
        """Return the 3D image mask array."""
        return self.mask_img.get_fdata()

    @property
    def img_arr(self):
        """Return the 3D image array."""
        return self.img.get_fdata()

    @property
    def resolution_xyz(self):
        """Return the delta x,y,z of the image."""
        return self.img.header["pixdim"][[1, 2, 3]]

    def map_coordinates(self, coord: Tuple, coord_type: str):
        """Save masked image to disc."""
        affine = self.img.affine

        # assign the affine transformation to use depending
        # if going from vox -> mm, or mm -> vox
        if coord_type == "vox":
            # from mm -> voxel
            inv_affine = npl.inv(affine)
        else:
            # from voxel -> mm
            inv_affine = affine

        new_coord = apply_affine(inv_affine, coord)
        return new_coord

    def apply_mask_xyz(self, coords_mm: Dict) -> Dict:
        """
        Filter out electrodes that do not fall within brain matter of a CT image.

        Parameters
        ----------
           elec_coords_mm: dict(str: ndarray)
               Dictionary that maps contact label to corresponding
               xyz coordinates.

        Returns
        -------
           elec_coords_mm: dict(str: ndarray)
               A dictionary of contact coordinates in CT voxels that fall
               within the binarized brain mask. The keys are individual
               contact labels, the values are the corresponding xyz coordinates
               in CT space.
        """
        # brain masked data has 0 outside the brain
        brainmasked_ct_img = self.get_masked_img()
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        brainmasked_ct_affine = brainmasked_ct_img.affine

        # Obtain inverse affine matrix to transform from xyz to CT voxel
        inv_affine = npl.inv(brainmasked_ct_affine)

        # Filter out channels not within mask
        _elec_coords_mm = {}
        for label, ch_coord in coords_mm.items():
            vox = apply_affine(inv_affine, ch_coord)
            idx = tuple(map(int, vox))
            if brainmasked_ct_data[idx] != 0:
                _elec_coords_mm[label] = ch_coord
        coords_mm = _elec_coords_mm

        return coords_mm

    def get_masked_img(self) -> nb.Nifti2Image:
        """
        Apply a binarized brain mask to an input CT image file.

        Returns
        -------
            masked_brain_img: NiBabel Image Object
                masked CT NiBabel Image Object in the same voxel space
                as input.
        """
        # obtain image arrays
        brain_img_arr = self.img_arr
        mask_arr = self.mask_arr

        # binarize mask to 0 and 1
        mask_arr[mask_arr > 0] = 1

        # multiply element wise brain array and mask array
        masked_brain = np.multiply(brain_img_arr, mask_arr)

        # return nibabel image
        masked_brain_img = nb.Nifti2Image(masked_brain, self.img.affine)
        return masked_brain_img

    def save_masked_img(self, fpath):
        """Save masked image to disc."""
        masked_brain_img = self.get_masked_img()
        masked_brain_img.to_filename(fpath)


class ClusteredBrainImage(BrainImage):
    """
    Brain image segmented into clusters.

    Inherits from :func:BrainImage
    """

    def __init__(self, *args):
        super(ClusteredBrainImage, self).__init__(*args)

    def _interpolate_points(self, clusters, num_to_fill, numcontacts):
        """
        Interpolate specified number of points.

        To fill between each cluster
        for a given electrode.

        Parameters
        ----------
            clusters: dict(str: ndarray)
                Dictionary of clustered_voxels for a given electrode.

            num_to_fill: ndarray
                Numpy array specifying the number of clustered_voxels to interpolate
                between each cluster centroid, i.e. num_to_fill[i] is the
                number of points to interpolate between the centroid of
                clustered_voxels[i] and clustered_voxels[i+1].

        Returns
        -------
            clustered_voxels: dict(str: ndarray)
                Updated dictionary of clustered_voxels for a given electrode, now
                containing interpolated points.
        """
        max_cluster_id = max(clusters.keys())

        # Obtain centroids
        centroids = self._compute_centroids(clusters)
        centroid_coords = np.array(list(centroids.values()))

        # Compute distance from each contact to its next immediate neighbor
        dists = np.diff(centroid_coords, axis=0)

        # Obtain dictionary of interpolated points
        for idx, num_to_interp in enumerate(num_to_fill):
            centroid = centroid_coords[idx]
            dist_to_next_centroid = dists[idx]
            new_pts = self._add_point(centroid, dist_to_next_centroid, num_to_interp)
            for new_pt in new_pts:
                # do not interpolate more clusters then number of contacts there are
                if len(clusters.keys()) >= numcontacts:
                    return clusters

                clusters[max_cluster_id + 1] = []
                clusters[max_cluster_id + 1].append(new_pt)
                max_cluster_id += 1

        return clusters

    def _add_point(self, centroid, dist_to_next_centroid, num_to_interp):
        """Add point in linear 3D space."""
        if centroid.shape != dist_to_next_centroid.shape:
            raise RuntimeError(
                "The centroid and distance should be " "in the same space."
            )

        # create linearly spaced proportion
        l1 = np.linspace(0, 1, num_to_interp)

        # add those proportions to the current centroid
        new_pts = centroid + dist_to_next_centroid * l1[:, None]

        return new_pts

    def _fill(self, centroids, num_to_fill, dists, max_cluster_id):
        """
        Insert new points between contacts.

        That are far apart for a given
        electrode and specified number of points to fill.

        Parameters
        ----------
            centroids: ndarray
                Numpy array of centroid coordinates along an electrode that
                need to have points interpolated between them. The centroids
                are ordered based on proximity to the first contact along the
                electrode, with the entry at index 0 being the contact closest
                to first contact.

            num_to_fill: ndarray
                Numpy array consisting of the number of points to interpolate
                between each centroid, i.e. num_to_fill[i] specifies the number
                of points to interpolate between centroids[i] and its immediate
                neighbor that is farther from the first contact along the
                electrode.

            dists: ndarray
                Numpy array consisting of the distances between each contact
                and its neighbor in the original spatial ordering, i.e.
                dists[i] is the L2 distance between centroids[i] and its
                immediate neighbor that is farther from the first contact along
                the electrode.

        Returns
        -------
            new_pts: dict(int: ndarray)
                Dictionary of new interpolated centroids to add.
        """
        new_pts = {}
        for i, num in enumerate(num_to_fill):
            for portion in range(1, num + 1):
                frac = portion / (num + 1)

                # Add fraction of the distance along direction to next
                # immediate neighbor
                new_pt = centroids[i] + (frac * dists[i])
                new_pts[max_cluster_id + portion] = np.array([new_pt])

                # Update the max_cluster_id to avoid overwriting any clustered_voxels
                max_cluster_id += num

        return new_pts

    def _identify_skull_voxel_clusters(
        self, voxel_clusters: Dict, skull_cluster_size: int = 200
    ):
        """
        Classify the abnormal clustered_voxels that are extremely large.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall.

                The keys of the dictionary are channel names,
                the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

        Returns
        -------
            skull_cluster_ids: dict(str: List[int])
                Dictionary of cluster ids thought to be large due to close
                proximity to the skull.
        """
        skull_clusters = []

        for cluster_id, voxels in voxel_clusters.items():
            # Average size of normal cluster is around 20-25 points
            cluster_size = len(voxels)

            # Skull clustered_voxels are likely to be very large
            if cluster_size > skull_cluster_size:
                skull_clusters.append(cluster_id)

        return skull_clusters

    def _identify_merged_voxel_clusters(self, voxel_clusters: Dict):
        """
        Classify the abnormal clustered_voxels that are extremely large.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

        Returns
        -------
            merged_cluster_ids: dict(str: List[int])
                Dictionary of cluster ids thought to be large due to lack of
                sufficient separation between two channels in image.
        """
        merged_cluster_ids = []

        for cluster_id, points in voxel_clusters.items():
            # Average size of normal cluster is around 20-25 points
            cluster_size = len(points)

            # Merged clustered_voxels are likely to be moderately large
            if 50 <= cluster_size <= 200:
                merged_cluster_ids.append(cluster_id)

        return merged_cluster_ids

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

    def _pare_cluster(self, points_in_cluster, qtile, centroid=None):
        """Pare a cluster down by a quantile around a centroid."""
        if centroid is None:
            # get the centroid of that cluster
            centroid = np.mean(points_in_cluster, keepdims=True)

        # compute spatial variance of the points inside cluster
        var = np.var(points_in_cluster, axis=0)

        # store the pared clsuters
        pared_cluster = []

        # Include points that have a z-score within specified quantile
        for pt in points_in_cluster:
            # Assuming the spatial distribution of points is
            # approximately Gaussian, the outermost channel will be
            # approximately the centroid of this cluster.
            diff = pt - centroid
            z = npl.norm(np.divide(diff, np.sqrt(var)))
            if norm.cdf(z) <= qtile:
                pared_cluster.append(pt)
        return pared_cluster

    def _pare_clusters_on_electrode(self, voxel_clusters, skull_cluster_ids, qtile):
        """
        Pare down skull clustered_voxels.

        Only considering points close to the
        centroid of the oversized cluster.

        Parameters
        ----------
            voxel_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

            skull_cluster_ids: dict(str: List[int])
                Dictionary of cluster ids thought to be large due to close
                proximity to the skull.

            qtile: float
                The upper bound quantile distance that we will consider for
                being a part of the resulting pared clustered_voxels.

        Returns
        -------
            voxel_clusters: dict(str: dict(str: ndarray))
                Dictionary of skull clustered_voxels that have been resized.
        """
        # for elec in skull_cluster_ids:
        # Get the coordinate for the outermost labeled channel from user in
        # last_chan_coord = list(sparse_labeled_contacts.values())[-1]
        for cluster_id in skull_cluster_ids:
            # get the clustered points for this cluster ID
            points_in_cluster = voxel_clusters[cluster_id]

            # pare the cluster based on quantile near the centroid
            pared_cluster = self._pare_cluster(points_in_cluster, qtile, centroid=None)

            # Sanity check that we still have a non-empty list
            if pared_cluster != []:
                voxel_clusters[cluster_id] = np.array(pared_cluster)

        return voxel_clusters

    def _compute_centroids(self, voxel_clusters):
        """
        Return the centroids of each channel label.

        Given a list of voxels per channel label.

        Parameters
        ----------
            chanxyzvoxels: dict()
                dictionary of electrodes and corresponding dictionary of
                channels to centroid xyz coordinates.

        Returns
        -------
            centroids: dict()
                dictionary of channels and corresponding centroid coordinates.
        """
        centroids = {}
        for label, voxels in voxel_clusters.items():
            centroids[label] = np.mean(np.array(voxels), axis=0)
        return centroids

    def _unfuse_clusters_on_entry_and_exit(
        self,
        voxel_clusters: Dict,
        merged_cluster_ids: List[int],
        contact_nums,
    ):
        """
        Unfuse merged clustered_voxels.

        By using KMeans clustering to split the large
        cluster into two.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are channel names,
                the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

            merged_cluster_ids: List[int]
                Dictionary of cluster ids thought to be large due to lack of
                sufficient separation between two channels in image.

        Returns
        -------
            voxel_clusters: dict(str: ndarray)
                Dictionary of skull clustered_voxels that have been resized.

        """
        # We need to be sure not to overwrite any clustered_voxels, so using values
        # greater than the max cluster ID will guarantee that
        max_cluster_id = max(voxel_clusters.keys())
        print(max_cluster_id)

        for cluster_id in merged_cluster_ids:
            # Use KMeans to separate cluster into two clustered_voxels
            cluster = voxel_clusters[cluster_id]
            kmeans = KMeans(n_clusters=2, random_state=0).fit(cluster)
            y = kmeans.labels_

            # Obtain points in each cluster
            cluster0, cluster1 = cluster[y == 0], cluster[y == 1]

            # Update the dictionary to separate the two clustered_voxels
            voxel_clusters[max_cluster_id + 1] = cluster0
            voxel_clusters[max_cluster_id + 2] = cluster1
            del voxel_clusters[cluster_id]

            max_cluster_id += 2

        return voxel_clusters

    def _order_clusters(self, clusters, first_contact_coord):
        """
        Order a dictionary of clustered_voxels.

        Based on distance of the centroid of
        the cluster from the contact labeled with the number 1 given from user
        input.

        Parameters
        ----------
            clusters: dict(str: ndarray)
                Dictionary of clustered_voxels for an electrode.

            first_contact: ndarray
                1x3 Numpy array of coordinates for the most medial labeled
                contact from user input.

        Returns
        -------
            sorted_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels for an electrode sorted based on
                proximity to the specified first_contact.
        """
        # compute centroids
        centroids = self._compute_centroids(clusters)

        # Sort based on proximity to most medial (innermost) contact
        sorted_centroids = sorted(
            centroids.items(), key=lambda x: npl.norm(x[1] - first_contact_coord)
        )
        sorted_centroids = dict(sorted_centroids)

        # Restore clustered_voxels, now in sorted order
        sorted_clusters = {k: clusters[k] for k in sorted_centroids}

        return sorted_clusters

    def compute_clusters_with_threshold(self, threshold: float = 0.630):
        """
        Apply a thresholding based algorithm and then running connected clustering using skimage.

        This will then return clustered_voxels of voxels that belong to distinct contact groups.
        http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label

        Parameters
        ----------
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
        brainmasked_CT_arr = self.get_masked_img().get_fdata()

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

        self.clusters = clusters

        return clusters, numobj

    def compute_cylindrical_clusters(
        self, clustered_voxels, entry_point_vox, exit_point_vox, radius
    ):
        """Compute clusters in cylinders."""
        # each electrode corresponds now to a separate cylinder
        voxel_clusters_in_cylinder = {}

        # go through each electrode apply cylinder filter based on entry/exit point
        # remove all voxels in clusters that are not within cylinder
        for cluster_id, cluster_voxels in clustered_voxels.items():
            # Obtain all points within the electrode endpoints pt1, pt2
            contained_voxels = []
            for voxel in cluster_voxels:
                if self._is_point_in_cylinder(
                    entry_point_vox, exit_point_vox, radius, voxel
                ):
                    contained_voxels.append(voxel)
            if not contained_voxels:
                continue

            # store all voxel clusters that are within this cylinder
            voxel_clusters_in_cylinder[cluster_id] = np.array(cluster_voxels)
        return voxel_clusters_in_cylinder

    def assign_sequential_labels(
        self, voxel_clusters, first_contact, first_contact_coord
    ):
        """
        Assign stereo-EEG electrode labels to clustered_voxels.

        Using labeling given by user.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels grouped by electrode.

            name: dict(str: ndarray)
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

        Returns
        -------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels grouped by electrode labeled with
                stereo-EEG labeling convention.
        """
        elec, ch_num = re.match("^([A-Za-z]+[']?)([0-9]+)$", first_contact).groups()
        ch_num = int(ch_num)

        # Ensure that the clustered_voxels are well-ordered
        # voxel_clusters = self._order_clusters(voxel_clusters, first_contact_coord)

        # Convert to numpy array to make indexing easier
        clusters = np.array(list(voxel_clusters.values()))

        if "'" in first_contact:
            # Electrode is inserted into left of brain, so we add in
            # forward order
            labeled_clusters = {
                elec + str(i + ch_num): clusters[i] for i in range(len(clusters))
            }
        else:
            # Electrode is inserted into left of brain, so we add in
            # reverse order
            labeled_clusters = {
                elec + str(i + ch_num): clusters[i]
                for i in range(len(clusters) - 1, -1, -1)
            }
        return labeled_clusters

    def fill_clusters_with_spacing(
        self,
        voxel_clusters: Dict,
        entry_ch: Contact,
        num_contacts,
        contact_spacing_mm: float,
    ):
        """
        Assist in filling in gaps in clustering.

        Compute the distances between a given channel and its immediate neighbor on the right.
        The last channel has a distance set to 0.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary with keys being electrode names and values being
                dictionaries consisting of entries of channel names and their
                corresponding centroid coordinates.

            gap_tolerance: float
                Max L2 distance (in voxel space) between adjacent centroids.

        Returns
        -------
            processed_clusters: dict(str: ndarray)
                Updated versions of final_xyz_centroids with adjusted labeling and
                dists with updated distances.
        """
        # Obtain labeled contacts given by user
        first_label = entry_ch.name
        first_coord = entry_ch.coord

        # Order clustered_voxels for an electrode based on proximity to the
        # innermost contact given by the user
        voxel_clusters = self._order_clusters(voxel_clusters, first_coord)

        # In rare case when there is only one cluster that belongs to a
        # given electrode, we exclusively use the user specified contacts
        if len(voxel_clusters) < 2:
            raise RuntimeError(
                "There is only one voxel cluster for " f"{entry_ch.name}?"
            )

        # Obtain centroid coordinates
        voxel_centroids = self._compute_centroids(voxel_clusters)

        # get array of the actual centroid coordinates
        centroid_coords = np.array(list(voxel_centroids.values()))

        # convert to xyz
        centroid_coords_mm = [
            self.map_coordinates(coord, coord_type="mm") for coord in centroid_coords
        ]
        # Compute L2 distance between adjacent centroids
        contact_distances = np.diff(centroid_coords_mm, axis=0)
        abs_distances = npl.norm(contact_distances, axis=1)

        # compute the median distance
        med_dist = np.median(abs_distances)
        contact_spacing_mm += max((med_dist - contact_spacing_mm) + 0.25, 0.5)

        # Compute the number of points to interpolate based on
        # specified gap tolerance
        num_to_fill = np.array(
            np.floor(abs_distances / contact_spacing_mm), dtype=np.uint16
        )
        print(f"Number to fill for {entry_ch}", num_to_fill)

        # Update the dictionary to include the new interpolated points
        voxel_clusters = self._interpolate_points(
            voxel_clusters, num_to_fill, len(num_contacts)
        )

        # Re-sort after adding interpolated points
        voxel_clusters = self._order_clusters(voxel_clusters, first_coord)

        return voxel_clusters

    def fill_gaps(self, voxel_clusters: Dict, gap_tolerance: float, entry_ch, exit_ch):
        """
        Assist in filling in gaps in clustering.

        Compute the distances between a given channel and its immediate neighbor on the right.
        The last channel has a distance set to 0.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary with keys being electrode names and values being
                dictionaries consisting of entries of channel names and their
                corresponding centroid coordinates.

            gap_tolerance: float
                Max L2 distance (in voxel space) between adjacent centroids.

        Returns
        -------
            processed_clusters: dict(str: ndarray)
                Updated versions of final_xyz_centroids with adjusted labeling and
                dists with updated distances.
        """
        # Obtain labeled contacts given by user
        first_label, first_coord = entry_ch
        last_label, last_coord = exit_ch

        # Order clustered_voxels for an electrode based on proximity to the
        # innermost contact given by the user
        voxel_clusters = self._order_clusters(voxel_clusters, first_coord)

        # In rare case when there is only one cluster that belongs to a
        # given electrode, we exclusively use the user specified contacts
        if len(voxel_clusters) < 2:
            raise RuntimeError(
                "There is only one voxel cluster for "
                f"{first_label} and {last_label}?"
            )
        # Obtain centroid coordinates
        voxel_centroids = self._compute_centroids(voxel_clusters)

        # get array of the actual centroid coordinates
        centroid_coords = np.array(list(voxel_centroids.values()))

        # Compute L2 distance between adjacent centroids
        dists = npl.norm(np.diff(centroid_coords, axis=0), axis=1)

        # Compute the number of points to interpolate based on
        # specified gap tolerance
        num_to_fill = np.floor(np.array(dists / gap_tolerance))

        print(centroid_coords)
        print(dists)
        print(f"Number to fill for {first_label} and {last_label}", num_to_fill)

        # Update the dictionary to include the new interpolated points
        voxel_clusters = self._interpolate_points(voxel_clusters, num_to_fill)

        # Re-sort after adding interpolated points
        voxel_clusters = self._order_clusters(voxel_clusters, first_coord)

        return voxel_clusters

    def bruteforce_correction(self, voxel_clusters, entry_ch, exit_ch):
        """
        Check for egregiously poor output by the algorithm.

        If it is large,
        use the brute force approach - use the provided contacts and
        interpolate specified number of evenly spaced poitns between the
        provided contacts.

        Parameters
        ----------
            voxel_clusters: dict(str: ndarray)
                Dictionary of clustered_voxels grouped by electrode labeled with
                stereo-EEG labeling convention.

        Returns
        -------
            labeled_clusters: dict(str: ndarray)
                Updated dictionary of clustered_voxels grouped by electrode labeled
                with stereo-EEG labeling convention with high error electrodes
                being replaced with the brute force approximation.
        """
        # Obtain labeled contacts given by user
        first_label, first_coord = entry_ch
        last_label, last_coord = exit_ch

        # Compute the number of points to linearly interpolate
        first_num = int(re.findall(r"\d+", first_label)[0])
        last_num = int(re.findall(r"\d+", last_label)[0])
        diff = last_num - first_num
        num_to_fill = np.array([diff - 1])

        # Compute centroids for this electrode
        centroids = self._compute_centroids(voxel_clusters)

        # Use the first and last centroid to compute error metric
        first_centroid = centroids.get(first_label, [])
        last_centroid = centroids.get(last_label, [])

        # Compute L2 distance between algorithm's prediction and provided
        # coordinates
        err_first = 0
        err_last = 0
        count = 0
        if len(first_centroid) > 0:
            err_first = npl.norm(first_centroid - first_coord)
            count += 1
        if len(last_centroid) > 0:
            err_last = npl.norm(last_centroid - last_coord)
            count += 1

        # Compute average of errors successfully computed
        if count == 0:
            err = 0
        else:
            err = (err_first + err_last) / count

        # Ignore errors smaller than 5.0
        if err <= 5.0:
            return voxel_clusters

        # Overwrite the algorithm's output for this electrode with the
        # user provided input
        voxel_clusters = {
            1: np.array([first_coord]),
            2: np.array([last_coord]),
        }

        # Interpolate the points and update the dictionary
        voxel_clusters = self._interpolate_points(voxel_clusters, num_to_fill)

        return voxel_clusters

    def cluster_2_centroids(self, clusters):
        """
        Compute and store the centroid for each cluster.

        Parameters
        ----------
            clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels grouped by electrode.
            final_xyz_centroids: dict()
                a dictionary with keys being electrode names and values being
                dictionaries consisting of entries of channel names and their
                corresponding centroid coordinates.
            dists: dict()
                a dictionary with keys being electrode names and values being
                dictionaries consisting of entries of channel names and their
                corresponding distance to their most immediate neighbor on their
                right. The rightmost channel for a given electrode has a distance
                set to 0.
            gap_tolerance: dict()
                maximum distance we will allow two adjacent contacts to be before
                no longer considering them adjacent.

        Returns
        -------
            centroids: dict(str: dict(str: ndarray))
                Dictionary of centroids grouped by electrode.
        """
        centroids = {}
        for elec in clusters:
            centroids[elec] = self._compute_centroids(clusters[elec])
        return centroids

    def vox_2_xyz(self, centroids_vox, affine):
        """
        Convert finalized dictionary of centroids from CT voxel coordinates to xyz coordinates.

        Parameters
        ----------
        centroids_vox: dict(str: dict(str: ndarray))
            Properly labeled dictionary of centroids in CT voxel
            coordinates.

        affine: ndarray
            Transformation matrix for applying affine transformation to
            coordinates to convert CT voxel coordinates to xyz coordinates.

        Returns
        -------
            centroids_xyz: dict(str: dict(str: ndarray))
                Properly labeled dictionary of centroids in xyz coordinates.
        """
        centroids_xyz = {}
        for elec in centroids_vox.keys():
            centroids_xyz[elec] = {}
            for chan in centroids_vox[elec]:
                centroids_xyz[elec][chan] = apply_affine(
                    affine, centroids_vox[elec][chan]
                )
        return centroids_xyz

    def _compute_voxels_around_centroid(
        self, voxel_centroid, voxel_clusters, memo_clusterids
    ):
        """Get voxels around a specific centroid within a ball-quantile."""
        # compute voxel cloud around this centroid
        new_cluster_voxels = []
        for _cluster_id in voxel_clusters.keys():
            if _cluster_id in memo_clusterids.keys():
                continue
            _cluster_voxels = voxel_clusters[_cluster_id]

            # is our new desired centroid within the convex hull?
            if _in_hull(_cluster_voxels, voxel_centroid):
                new_cluster_voxels = self._pare_cluster(
                    _cluster_voxels, qtile=0.5, centroid=voxel_centroid
                )
        return new_cluster_voxels

    def _search_clusters_for_point(
        self, this_elec_xyz, approx_new_centroid, lower_bound_point, upper_bound_point
    ):
        """Search a set of clusters for a 3D point."""
        matched_ids = []

        centroid_dists = []
        for cluster_id, point_cluster in this_elec_xyz.items():
            # if _in_hull(point_cluster, approx_new_centroid):
            #     matched_ids.append(cluster_id)
            centroid_dists.append(np.mean(point_cluster, axis=0) - approx_new_centroid)
            matched_ids.append(cluster_id)

        print(np.array(centroid_dists).shape)
        abs_distances = np.linalg.norm(centroid_dists, axis=0)
        print(abs_distances.shape)
        matched_ind = np.argmin(abs_distances)
        matched_ids = [matched_ids[matched_ind]]

        if len(matched_ids) > 1:
            print(
                "Multiple clusters satisfying this condition! "
                "Returning the first one found."
            )
        if matched_ids == []:
            return None

        return this_elec_xyz[matched_ids[0]]

    def correct_labeled_clusters(
        self,
        this_elec_voxels: Dict,
        entry_ch: Contact,
        exit_ch: Contact,
        contact_spacing_mm: float,
    ):
        """
        Recursive iteration of electrodes and correction.

        The basic idea is that we go through each cluster set,

        Parameters
        ----------
        this_elec_voxels :
        entry_ch :
        exit_ch :
        contact_spacing_mm :
        """
        # allow slack of the contact spacing
        epsilon = 0.5

        # generate unit vector from entry to exit channel
        unit_vec = np.array(
            _compute_vector(entry_ch.coord, exit_ch.coord, unitSphere=True)
        )

        # go through each cluster set of points
        current_centroid = entry_ch.coord

        # memoize the cluster ids
        memo_clusterids = dict()

        for cluster_id, cluster_voxels in this_elec_voxels.items():
            # compute 3D voxel centroid
            voxel_centroid = np.mean(cluster_voxels, axis=0)

            assert len(voxel_centroid) == 3

            # compute distance between centroids
            distance = np.linalg.norm(voxel_centroid - current_centroid)

            # outside the epsilon-ball allowed for contact spacing
            if (distance > contact_spacing_mm + epsilon) or (
                distance < contact_spacing_mm - epsilon
            ):
                # compute new centroid
                voxel_centroid = current_centroid + contact_spacing_mm * unit_vec

                # compute a new set of clustered voxels based on centroid desired
                new_cluster_voxels = self._compute_voxels_around_centroid(
                    voxel_centroid, this_elec_voxels, memo_clusterids
                )

                # if empty, manually set one centroid
                if new_cluster_voxels == []:
                    new_cluster_voxels = [voxel_centroid]

                # set the newly identified cluster voxels
                this_elec_voxels[cluster_id] = new_cluster_voxels

            # update the current centroid
            current_centroid = voxel_centroid

            # update memoized clusters we've checked
            memo_clusterids[cluster_id] = 1

        # on the exit contact, check that the manually labeled centroid is within the cluster
        if not _in_hull(this_elec_voxels[cluster_id], exit_ch.coord):
            this_elec_voxels[cluster_id] = [exit_ch.coord]

        return this_elec_voxels

    def bruteforce_correctionv2(
        self,
        this_elec_xyz: Dict,
        entry_ch: Contact,
        exit_ch: Contact,
        contact_spacing_mm: float,
        num_contacts,
    ):
        """
        Recursive iteration of electrodes and correction.

        The basic idea is that we go through each cluster set,

        Parameters
        ----------
        this_elec_xyz :
        entry_ch :
        exit_ch :
        contact_spacing_mm :
        """
        # allow slack of the contact spacing
        epsilon = 0.5

        # generate unit vector from entry to exit channel
        unit_vec = np.array(
            _compute_vector(entry_ch.coord, exit_ch.coord, unitSphere=True)
        )

        # create a copy of the dictionary set of cluster points
        _cluster_xyz = this_elec_xyz.copy()
        cluster_ids = list(this_elec_xyz.keys())
        # create new copy of electrode cluster points
        new_cluster_xyz = collections.defaultdict(list)

        # loop until we've filled a point w/ each contact
        prev_centroid = entry_ch.coord
        i = 0
        while i < num_contacts:
            lower_bound_point = prev_centroid + unit_vec * (
                contact_spacing_mm - epsilon
            )
            upper_bound_point = prev_centroid + unit_vec * (
                contact_spacing_mm + epsilon
            )
            approx_new_centroid = prev_centroid + unit_vec * contact_spacing_mm
            cluster_points = self._search_clusters_for_point(
                this_elec_xyz, approx_new_centroid, lower_bound_point, upper_bound_point
            )
            if cluster_points is None:
                prev_centroid = approx_new_centroid
                cluster_points = [prev_centroid]
            else:
                prev_centroid = np.mean(cluster_points, axis=0)

            # create new set of clusters
            new_cluster_xyz[i] = cluster_points
            i += 1
        return new_cluster_xyz

        # # perform loop indefinitely until constraints are met
        # while len(this_elec_xyz.keys()) <= num_contacts and \
        #         len(memo_clusterids.keys()) != len(this_elec_xyz.keys()):
        #     cluster_id = cluster_ids[i]
        #     # go through the identified clusters so far
        #     cluster_voxels = this_elec_xyz[cluster_id]
        #     # compute 3D voxel centroid
        #     voxel_centroid = np.mean(cluster_voxels, axis=0)
        #     # compute distance between centroids
        #     distance = np.linalg.norm(voxel_centroid - current_centroid)
        #
        #     # outside the epsilon-ball allowed for contact spacing
        #     if distance < epsilon:
        #         # compute new centroid
        #         voxel_centroid = current_centroid + contact_spacing_mm * unit_vec
        #
        #         # compute a new set of clustered voxels based on centroid desired
        #         new_cluster_voxels = self._compute_voxels_around_centroid(voxel_centroid,
        #                                                                   this_elec_xyz,
        #                                                                   memo_clusterids)
        #         # if empty, manually set one centroid
        #         if new_cluster_voxels == []:
        #             new_cluster_voxels = [voxel_centroid]
        #
        #         # set the newly identified cluster voxels
        #         new_cluster_xyz[cluster_id] = new_cluster_voxels
        #     else:
        #         new_cluster_xyz[cluster_id] = cluster_voxels
        #
        #     # update the current centroid
        #     current_centroid = voxel_centroid
        #     # update memoized clusters we've checked
        #     memo_clusterids[cluster_id] = 1
        #     i+=1
        #
        # # on the exit contact, check that the manually labeled centroid is within the cluster
        # if not _in_hull(this_elec_xyz[cluster_id], exit_ch.coord):
        #     this_elec_xyz[cluster_id] = [exit_ch.coord]
        #
        # return this_elec_xyz
