import collections
import re
from typing import Dict, List

import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine
from scipy.stats import norm
from sklearn.cluster import KMeans


class PostProcessor:
    """Class of postprocessor grouping functions."""

    @classmethod
    def _identify_skull_voxel_clusters(
        cls, electrode_voxel_clusters: Dict, skull_cluster_size: int = 200
    ):
        """
        Classify the abnormal clustered_voxels that are extremely large.

        Parameters
        ----------
            electrode_voxel_clusters: dict(str: ndarray)
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

        for cluster_id, voxels in electrode_voxel_clusters.items():
            # Average size of normal cluster is around 20-25 points
            cluster_size = len(voxels)

            # Skull clustered_voxels are likely to be very large
            if cluster_size > skull_cluster_size:
                skull_clusters.append(cluster_id)

        return skull_clusters

    @classmethod
    def _identify_merged_voxel_clusters(self, electrode_voxel_clusters: Dict):
        """
        Classify the abnormal clustered_voxels that are extremely large.

        Parameters
        ----------
            electrode_voxel_clusters: dict(str: ndarray)
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

        for cluster_id, points in electrode_voxel_clusters.items():
            # Average size of normal cluster is around 20-25 points
            cluster_size = len(points)

            # Merged clustered_voxels are likely to be moderately large
            if 50 <= cluster_size <= 200:
                merged_cluster_ids.append(cluster_id)

        return merged_cluster_ids

    @classmethod
    def _pare_clusters_on_electrode(self, electrode_clusters, skull_cluster_ids, qtile):
        """
        Pare down skull clustered_voxels.

        Only considering points close to the
        centroid of the oversized cluster.

        Parameters
        ----------
            electrode_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

            skull_cluster_ids: dict(str: List[int])
                Dictionary of cluster ids thought to be large due to close
                proximity to the skull.

            sparse_labeled_contacts: dict(str: dict(str: ndarray))
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

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
            points_in_cluster = electrode_clusters[cluster_id]

            # get the centroid of that cluster
            cluster_centroid = np.mean(points_in_cluster, keepdims=True)

            var = np.var(points_in_cluster, axis=0)

            # store the pared clsuters
            pared_cluster = []

            # Include points that have a z-score within specified quantile
            for pt in points_in_cluster:
                # Assuming the spatial distribution of points is
                # approximately Gaussian, the outermost channel will be
                # approximately the centroid of this cluster.
                diff = pt - cluster_centroid
                z = npl.norm(np.divide(diff, np.sqrt(var)))
                if norm.cdf(z) <= qtile:
                    pared_cluster.append(pt)

            # Sanity check that we still have a non-empty list
            if pared_cluster != []:
                electrode_clusters[cluster_id] = np.array(pared_cluster)

        return electrode_clusters

    @classmethod
    def _unfuse_clusters_on_electrode(
        self, electrode_clusters: Dict, merged_cluster_ids: List[int]
    ):
        """
        Unfuse merged clustered_voxels.

        By using KMeans clustering to split the large
        cluster into two.

        Parameters
        ----------
            electrode_clusters: dict(str: ndarray)
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
        max_cluster_id = max(electrode_clusters.keys())

        for cluster_id in merged_cluster_ids:
            # Use KMeans to separate cluster into two clustered_voxels
            cluster = electrode_clusters[cluster_id]
            kmeans = KMeans(n_clusters=2, random_state=0).fit(cluster)
            y = kmeans.labels_

            # Obtain points in each cluster
            cluster0, cluster1 = cluster[y == 0], cluster[y == 1]

            # Update the dictionary to separate the two clustered_voxels
            electrode_clusters[max_cluster_id + 1] = cluster0
            electrode_clusters[max_cluster_id + 2] = cluster1
            del electrode_clusters[cluster_id]

            max_cluster_id += 2

        return electrode_clusters

    @classmethod
    def process_abnormal_clusters(
        self, cyl_filtered_clusters, sparse_labeled_clusters, qtile=0.875
    ):
        """
        Truncate the clustered_voxels closest to the skull.

        Which tend to be oversized, and separates clustered_voxels that
        appear to have grouped two
        contacts together. The cluster is filtered with 90% quantile filtering,
        using a coordinate from user-input as the mean for the cluster.

        Parameters
        ----------
            cyl_filtered_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels sorted by the cylinder/electrode
                in which they fall. The keys of the dictionary are electrode
                labels, the values of the dictionary are the cluster points
                from the threshold-based clustering algorithm that fell
                into a cylinder.

            grouped_electrode_voxels: dict(str: dict(str: ndarray))
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

            qtile: float
                The upper bound quantile distance that we will consider for
                being a part of the resulting pared clustered_voxels.

        Returns
        -------
            voxel_clusters: dict(str: dict(str: ndarray))
                Updated dictionary of clustered_voxels found along each electrode
                with the skull clustered_voxels having been resized through quantile
                filtering and large clustered_voxels having been separated.

        """
        # Identify abnormally merged clustered_voxels
        merged_cluster_ids = self._identify_merged_voxel_clusters(cyl_filtered_clusters)

        # identify possibly oversize voxel clusters at the skull within each cylinder
        skull_cluster_ids = self._identify_skull_voxel_clusters(cyl_filtered_clusters)

        print(f"Oversized clustered_voxels: {skull_cluster_ids}")
        print(f"Merged clustered_voxels: {merged_cluster_ids}")

        # Resize skull clustered_voxels
        cyl_filtered_clusters = self._pare_clusters_on_electrode(
            cyl_filtered_clusters, skull_cluster_ids, sparse_labeled_clusters, qtile
        )

        # Separate fused clustered_voxels into two new clustered_voxels and store them
        # in a dictionary
        cyl_filtered_clusters = self._unfuse_clusters_on_electrode(
            cyl_filtered_clusters, merged_cluster_ids, sparse_labeled_clusters
        )

        return cyl_filtered_clusters

    @classmethod
    def _compute_centroids(self, chanxyzvoxels):
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
        for label, coords in chanxyzvoxels.items():
            centroids[label] = np.mean(np.array(coords), axis=0)
        return centroids

    @classmethod
    def _order_clusters(self, clusters, first_contact):
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
            centroids.items(), key=lambda x: npl.norm(x[1] - first_contact)
        )
        sorted_centroids = collections.OrderedDict(sorted_centroids)

        # Restore clustered_voxels, now in sorted order
        sorted_clusters = {k: clusters[k] for k in sorted_centroids}

        return sorted_clusters

    @classmethod
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

    @classmethod
    def _interpolate_points(self, clusters, num_to_fill):
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
        idxs = np.nonzero(num_to_fill)

        # Obtain centroids
        centroids = self._compute_centroids(clusters)
        centroid_coords = np.array(list(centroids.values()))

        # Compute distance from each contact to its next immediate neighbor
        dists = np.diff(centroid_coords, axis=0)

        # Obtain dictionary of interpolated points
        new_pts = self._fill(
            centroid_coords[idxs], num_to_fill[idxs], dists[idxs], max_cluster_id
        )

        # Merge the dictionary of interpolated points into clustered_voxels
        clusters.update(new_pts)

        return clusters

    @classmethod
    def fill_gaps(
        self, electrode_clusters: Dict, gap_tolerance: float, entry_ch, exit_ch
    ):
        """
        Assist in filling in gaps in clustering.

        Compute the distances between a given channel and its immediate neighbor on the right.
        The last channel has a distance set to 0.

        Parameters
        ----------
            electrode_clusters: dict(str: ndarray)
                Dictionary with keys being electrode names and values being
                dictionaries consisting of entries of channel names and their
                corresponding centroid coordinates.

            gap_tolerance: float
                Max L2 distance (in voxel space) between adjacent centroids.

            sparse_labeled_contacts: dict(str: dict(str: ndarray))
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

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
        electrode_clusters = self._order_clusters(electrode_clusters, first_coord)

        # In rare case when there is only one cluster that belongs to a
        # given electrode, we exclusively use the user specified contacts
        if len(electrode_clusters) < 2:
            max_id = max(electrode_clusters.keys())

            electrode_clusters = {
                max_id + 1: np.array([first_coord]),
                max_id + 2: np.array([last_coord]),
            }
            first_num = int(re.findall(r"\d+", first_label)[0])
            last_num = int(re.findall(r"\d+", last_label)[0])
            diff = last_num - first_num
            num_to_fill = np.array([diff - 1])

        else:
            # Obtain centroid coordinates
            centroids = self._compute_centroids(electrode_clusters)
            centroid_coords = np.array(list(centroids.values()))

            # Compute L2 distance between adjacent centroids
            dists = npl.norm(np.diff(centroid_coords, axis=0), axis=1)

            # Compute the number of points to interpolate based on
            # specified gap tolerance
            num_to_fill = np.array(dists // gap_tolerance, dtype=np.uint16)

        # Update the dictionary to include the new interpolated points
        electrode_clusters = self._interpolate_points(electrode_clusters, num_to_fill)

        # Re-sort after adding interpolated points
        electrode_clusters = self._order_clusters(electrode_clusters, first_coord)

        return electrode_clusters

    @classmethod
    def assign_labels(self, final_clusters, sparse_labeled_contacts):
        """
        Assign stereo-EEG electrode labels to clustered_voxels.

        Using labeling given by user.

        Parameters
        ----------
            final_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels grouped by electrode.

            sparse_labeled_contacts: dict(str: dict(str: ndarray))
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

        Returns
        -------
            voxel_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels grouped by electrode labeled with
                stereo-EEG labeling convention.
        """
        labeled_clusters = {}
        for elec in final_clusters:
            # Obtain innermost contact from user specified contacts
            chs_on_elec = sparse_labeled_contacts[elec]
            ch_names = list(chs_on_elec.keys())
            first_contact = sparse_labeled_contacts[elec][ch_names[0]]

            # Ensure that the clustered_voxels are well-ordered
            final_clusters[elec] = self._order_clusters(
                final_clusters[elec], first_contact
            )

            # Convert to numpy array to make indexing easier
            clusters = np.array(list(final_clusters[elec].values()))

            if "'" in elec:
                # Electrode is inserted into left of brain, so we add in
                # forward order
                labeled_clusters[elec] = {
                    elec + str(i + 1): clusters[i] for i in range(len(clusters))
                }
            else:
                # Electrode is inserted into left of brain, so we add in
                # reverse order
                labeled_clusters[elec] = {
                    elec + str(i + 1): clusters[i]
                    for i in range(len(clusters) - 1, -1, -1)
                }
        return labeled_clusters

    @classmethod
    def bruteforce_correction(self, labeled_clusters, sparse_labeled_contacts):
        """
        Check for egregiously poor output by the algorithm.

        If it is large,
        use the brute force approach - use the provided contacts and
        interpolate specified number of evenly spaced poitns between the
        provided contacts.

        Parameters
        ----------
            labeled_clusters: dict(str: dict(str: ndarray))
                Dictionary of clustered_voxels grouped by electrode labeled with
                stereo-EEG labeling convention.

            sparse_labeled_contacts: dict(str: dict(str: ndarray))
                Sparse dictionary of labeled channels from user input. This
                dictionary contains exactly two channels for each electrode.

        Returns
        -------
            voxel_clusters: dict(str: dict(str: ndarray))
                Updated dictionary of clustered_voxels grouped by electrode labeled
                with stereo-EEG labeling convention with high error electrodes
                being replaced with the brute force approximation.
        """
        for elec in sparse_labeled_contacts:
            # Obtain user provided contacts
            first, last = list(sparse_labeled_contacts[elec].items())
            first_label, first_coord = first
            last_label, last_coord = last

            # Compute centroids for this electrode
            centroids = self._compute_centroids(labeled_clusters[elec])

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
                continue

            # Overwrite the algorithm's output for this electrode with the
            # user provided inpute
            labeled_clusters[elec] = {
                1: np.array([first_coord]),
                2: np.array([last_coord]),
            }

            # Compute the number of points to linearly interpolate
            first_num = int(re.findall(r"\d+", first_label)[0])
            last_num = int(re.findall(r"\d+", last_label)[0])
            diff = last_num - first_num
            num_to_fill = np.array([diff - 1])

            # Interpolate the points and update the dictionary
            labeled_clusters[elec] = self._interpolate_points(
                labeled_clusters[elec], num_to_fill
            )

        # Reassign SEEG labels
        labeled_clusters = self.assign_labels(labeled_clusters, sparse_labeled_contacts)

        return labeled_clusters

    @classmethod
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

    @classmethod
    def vox_2_xyz(self, centroids_vox, affine):
        """
        Convert finalized dictionary of centroids from CT voxel coordinates to xyz coordinates.

        Parameters
        ----------
            centroids_vox: dict(str: dict(str: ndarray))
                Properly labeled dictioanry of centroids in CT voxel
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
