import copy

import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine
from scipy.stats import norm
from sklearn.decomposition import PCA


class PostProcessor:
    @classmethod
    def _group_by_cylinder(self, cylindrical_filter_clusters):
        """
        Groups clusters by electrodes according to cylindrical boundaries.

        Parameters
        ----------
        cylindrical_filter_clusters: dict()
            Dictionary of clusters that have been cylindrically filtered.
            All points outside the cylindrical filter are zeroed out.
        """
        grouped_clusters = {}
        for electrode in cylindrical_filter_clusters:
            grouped_clusters[electrode] = []
            for cluster_id in cylindrical_filter_clusters[electrode]:
                for i in range(len(cylindrical_filter_clusters[electrode][cluster_id])):
                    point = cylindrical_filter_clusters[electrode][cluster_id][i]
                    grouped_clusters[electrode].append(point)
        return grouped_clusters

    @classmethod
    def _label_cylinder_clusters(self, grouped_cylindrical_clusters, clusters):
        """
        Parameters
        ----------
        grouped_cylindrical_clusters: dict()
            Nested dictionary of contacts and their point clouds
            grouped by electrode.

        clusters: dict()
            Dictionary of cluster ids to point clouds.

        Returns
        -------
        Cluster labeling replaced with channel labels rather than
        numerical id's
        """
        cluster_labels_per_electrode = {}
        for electrode in grouped_cylindrical_clusters:
            cluster_labels_per_electrode[electrode] = []
            for cluster_id in clusters:
                # Convert every cluster into tuples to be hashable
                a = map(tuple, clusters[cluster_id])
                b = map(tuple, grouped_cylindrical_clusters[electrode])
                intxn = list(set(a) & set(b))
                if len(intxn) > 0:
                    cluster_labels_per_electrode[electrode].append(cluster_id)
        return cluster_labels_per_electrode

    @classmethod
    def _get_cluster_sizes(self, grouped_cylindrical_clusters, clusters):
        """
        """
        cluster_sizes_by_electrode = {}
        for electrode in grouped_cylindrical_clusters:
            cluster_sizes_by_electrode[electrode] = []
            for cluster_id in grouped_cylindrical_clusters[electrode]:
                cluster = clusters[cluster_id]
                cluster_sizes_by_electrode[electrode].append(len(cluster))
        return cluster_sizes_by_electrode

    @classmethod
    def _typeify_abnormalities(
        self, grouped_cylindrical_clusters, cluster_sizes_by_electrode
    ):
        merged_clusters, skull_clusters = [], []
        for electrode, cluster_sizes in cluster_sizes_by_electrode.items():
            for i in range(len(cluster_sizes)):
                # Hyperparameters: threshold for which we consider a cluster
                # to be a skull cluster or fused cluster
                if cluster_sizes[i] > 200:
                    skull_clusters.append(grouped_cylindrical_clusters[electrode][i])
                elif cluster_sizes[i] >= 50 and cluster_sizes[i] <= 200:
                    merged_clusters.append(grouped_cylindrical_clusters[electrode][i])
        return merged_clusters, skull_clusters

    @classmethod
    def process_abnormal_clusters(
        self,
        clusters,
        elec_in_brain,
        cylindrical_filter_clusters,
        sparse_elec_labels,
        qtile=0.875,
    ):
        """
        Truncates the clusters closest to the skull, which tend to be 
        oversized, and separates clusters that appear to have grouped two 
        contacts together. The cluster is filtered with 90% quantile filtering, 
        using a coordinate from user-input as the mean for the cluster.

        Parameters
        –---------
        clusters: dict()
            Dictionary of clusters found at a given threshold.

        elec_in_brain: dict()
            Dictionary of all the contacts found within the brain.

        cylindrical_filter_clusters: dict()
            Dictionary of points first sorted by the cylindrical electrode
            in which they fall and then the cluster in which they fall.
            The keys of the dictionary are electrode labels, the values of
            the dictionary are the cluster points from the threshold-based
            clustering algorithm that fell into a cylinder.

        sparse_elec_labels: dict()
            sparse dictionary of electrode labels given by user input

        Returns
        -------
            Updated dictionary of clusters found along each electrode
            with the skull clusters having been resized through quantile 
            filtering and large clusters having been separated.

        """
        # Group the cylindrically bounded clusters by electrodes
        grouped_cylindrical_clusters = self._group_by_cylinder(
            cylindrical_filter_clusters
        )
        grouped_cylindrical_clusters = self._label_cylinder_clusters(
            grouped_cylindrical_clusters, clusters
        )
        cluster_sizes_by_electrode = self._get_cluster_sizes(
            grouped_cylindrical_clusters, clusters
        )

        # Identify abnormal clusters
        merged_clusters, skull_clusters = self._typeify_abnormalities(
            grouped_cylindrical_clusters, cluster_sizes_by_electrode
        )
        result = copy.deepcopy(cylindrical_filter_clusters)

        print(f"Oversized clusters: {skull_clusters}")
        print(f"Merged clusters: {merged_clusters}")

        # Resize skull clusters
        resized_clusters = {}
        for cluster_id in skull_clusters:
            cluster_points = clusters[cluster_id]
            var = np.var(cluster_points, axis=0)
            for electrode in grouped_cylindrical_clusters:
                if cluster_id in grouped_cylindrical_clusters[electrode]:
                    last_contact_label = sparse_elec_labels[electrode][-1]
                    last_contact = elec_in_brain[last_contact_label]
                    resized_clusters[cluster_id] = []
                    for point in cluster_points:
                        z = npl.norm(np.divide((point - last_contact), np.sqrt(var)))
                        # Only accept up to 87.5% quantile
                        if norm.cdf(z) <= qtile:
                            resized_clusters[cluster_id].append(point)
                    if not resized_clusters[cluster_id]:
                        del resized_clusters[cluster_id]

        # Update resulting dictionary after resizing skull clusters
        for electrode in result:
            for cluster_id in resized_clusters:
                if cluster_id in result[electrode].keys():
                    result[electrode][cluster_id] = resized_clusters[cluster_id]

        # Separate fused clusters into two new clusters and store them
        # in a dictionary
        unmerged_clusters = {}
        max_cluster_id = max(clusters.keys())
        for cluster_id in merged_clusters:
            pts = clusters[cluster_id]
            var = np.var(pts, axis=0)
            var = npl.norm(var)
            for electrode in grouped_cylindrical_clusters:
                if cluster_id in grouped_cylindrical_clusters[electrode]:
                    if len(sparse_elec_labels[electrode]) < 2:
                        raise ValueError(
                            f"Unable to split merged cluster for "
                            f"electrode {electrode}. "
                            f"Please manually label "
                            f"the second innermost channel."
                        )
                    first_contact_label = sparse_elec_labels[electrode][0]
                    second_contact_label = sparse_elec_labels[electrode][1]
                    first_contact = elec_in_brain[first_contact_label]
                    second_contact = elec_in_brain[second_contact_label]
                    unmerged_clusters[electrode] = {}
                    unmerged_clusters[electrode][max_cluster_id + 1] = []
                    unmerged_clusters[electrode][max_cluster_id + 2] = []
                    for point in pts:
                        z1 = npl.norm(np.divide((point - first_contact), np.sqrt(var)))
                        z2 = npl.norm(np.divide((point - second_contact), np.sqrt(var)))
                        if norm.cdf(z1) <= qtile:
                            unmerged_clusters[electrode][max_cluster_id + 1].append(
                                point
                            )
                        if norm.cdf(z2) <= qtile:
                            unmerged_clusters[electrode][max_cluster_id + 2].append(
                                point
                            )
                    max_cluster_id += 2

        # Update cluster dictionary
        for electrode in result:
            # delete the merged cluster id's
            for cluster_id in merged_clusters:
                if cluster_id in result[electrode].keys():
                    del result[electrode][cluster_id]
            # add the corresponding unmerged cluster id's
            if electrode in unmerged_clusters.keys():
                for cluster_id in unmerged_clusters[electrode]:
                    result[electrode][cluster_id] = unmerged_clusters[electrode][
                        cluster_id
                    ]

        relabeled_clusters = {}
        for electrode in result:
            num = 1
            relabeled_clusters[electrode] = {}
            for cluster_id in result[electrode]:
                new_id = electrode + str(num)
                relabeled_clusters[electrode][new_id] = result[electrode][cluster_id]
                num += 1
        return relabeled_clusters

    @classmethod
    def fill_one_point(self, p1, p2):
        """
        Computes midpoint.

        Parameters
        –---------
        p1: np.ndarray of dimension 1 x 3
            first point to compute midpoint.

        p2: np.ndarray of dimension 1 x 3
            second point to compute midpoint.

        Returns
        -------
            midpoint of the two specified points.
        """
        return (p1 + p2) / 2

    @classmethod
    def fill_two_points(self, p1, p2):
        """
        Computes trisection points.

        Parameters
        ----------
        p1: np.ndarray of dimension 1 x 3
            first point to compute trisection points.

        p2: np.ndarray of dimension 1 x 3
            second point to compute trisection points.

        Returns
        -------
            Two points which trisect the line segment formed by p1 and p2.
        """
        return ((2 / 3) * p1 + (1 / 3) * p2, (1 / 3) * p1 + (2 / 3) * p2)

    @classmethod
    def fill_three_points(self, p1, p2):
        """
        Computes quadrisection points.

        Parameters
        –---------
        p1: np.ndarray of dimension 1 x 3
            first point to compute quadrisection points.

        p2: np.ndarray of dimension 1 x 3
            second point to compute qaudrisection points.

        Returns
        -------
            Three points which quadrisect the line segment formed by p1 and p2.
        """
        return (
            (1 / 4) * p1 + (3 / 4) * p2,
            self.fill_one_point(p1, p2),
            (3 / 4) * p1 + (1 / 4) * p2,
        )

    @classmethod
    def shift_downstream_labels(self, electrode_name, idx, shift_factor, chan_dict):
        """
        Helps to relabel downstream labels from a given index to maintain
        sorted order.

        Parameters
        ----------
        electrode_name: str
            string that contains electrode name (e.g. L')

        idx: int
            index from which shifting should start

        shift_factor: int
            how much to shift downstream labels by

        chan_dict: dict
            dictionary of channels with standard labeling and their
            corresponding coordinates.

        Returns
        -------
            Updated version of chan_dict with shifted labels to allow
            for easy insertion.
        """
        start = len(list(chan_dict.keys()))
        for i in range(start, idx - 1, -1):
            cur_id = electrode_name + str(i)
            new_id = electrode_name + str(i + shift_factor)
            chan_dict[new_id] = chan_dict[cur_id]
            del chan_dict[cur_id]
        return chan_dict

    @classmethod
    def compute_dists(self, centroid_dict):
        """
        Compute the distances between a given channel and its immediate
        neighbor on the right. The last channel has a distance set to 0.

        Parameters
        ----------
        centroid_dict: dict()
            a dictionary with keys being electrode names and values being
            dictionaries consisting of entries of channel names and their
            corresponding centroid coordinates.

        Returns
        -------
            a dictionary with keys being electrode names and values being
            dictionaries consisting of entries of channel names and their
            corresponding distances as described above.
        """

        dists = {elec: {} for elec in list(centroid_dict.keys())}
        for electrode in centroid_dict:
            channels = list(centroid_dict[electrode].keys())
            num_chans = len(channels)
            centrs = centroid_dict[electrode]
            for i in range(num_chans - 1):
                cur = np.array(centrs[channels[i]])
                nxt = np.array(centrs[channels[i + 1]])
                dists[electrode][channels[i]] = npl.norm(cur - nxt)
            dists[electrode][channels[-1]] = 0
        return dists

    @classmethod
    def _strip_nans(self, channels):
        """
        Strip any nan values from a channel array.

        Parameters
        ----------
        channels: np.ndarray of dimension 1 x n
            numpy array of channels and possibly nan values.

        Returns
        -------
            numpy array with the same channels excluding nan values.
        """
        channels = np.array(channels)
        temp = [chan for chan in channels if not np.isnan(chan).all()]
        return temp

    @classmethod
    def reassign_labels(self, centroid_dict):
        """
        Assigns labels which follow the standard labeling convention for SEEG
        electrodes.

        Parameters
        ----------
            centroid_dict: dict()
                A dictionary where the keys are each electrode name and values
                are dictionaries with entries of channels (does not necessarily
                have to follow standard labeling convention) and their
                corresponding coordinates.

        Returns
        -------
            A dictionary where channels are correctly assigned and sorted in
            order. If the electrode name has a ', then the first entry in the
            dictionary will be the smallest number (usually 1) and every
            subsequent entry will be larger. Otherwise, the first entry
            in the dictionary will be the largest number and every subsequent
            entry will be smaller.
        """
        pca = PCA()
        result = {}
        for electrode in centroid_dict:
            result[electrode] = {}
            chans = np.array(list(centroid_dict[electrode].values()))
            stripped_chans = self._strip_nans(chans)
            pca.fit(stripped_chans)
            centroids_pca = pca.transform(stripped_chans)
            centroids_new = pca.inverse_transform(centroids_pca)
            sorted_idxs = np.argsort(centroids_new[:, 0])
            sorted_orig_ids = np.array(list(centroid_dict[electrode].keys()))[
                sorted_idxs
            ]

            # assert len(sorted_idxs) == len(centroid_dict[electrode].keys())
            # Left side of the brain
            if electrode[-1] == "'":
                for i, chan in enumerate(sorted_idxs):
                    new_id = electrode + str(i + 1)
                    result[electrode][new_id] = centroid_dict[electrode][
                        sorted_orig_ids[i]
                    ]
            # Right side of the brain
            else:
                for i, chan in enumerate(sorted_idxs):
                    new_id = electrode + str(len(sorted_orig_ids) - i)
                    result[electrode][new_id] = centroid_dict[electrode][
                        sorted_orig_ids[i]
                    ]
        return result

    @classmethod
    def fill_gaps(self, final_centroids, dists, gap_tolerance):
        """
        Assist in filling in gaps in clustering.

        Parameters
        –---------
        final_centroids: dict()
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
            Updated versions of final_centroids with adjusted labeling and
            dists with updated distances.
        """
        for electrode in dists:
            for i, chan in enumerate(dists[electrode]):
                if i == len(list(dists[electrode].keys())) - 1:
                    continue
                chan_list = np.array(list(dists[electrode].keys()))
                # Hyperparameters: What distance to consider having skipped
                # one, two, or three channels
                # Need to fill in one channel
                if (
                    dists[electrode][chan] > gap_tolerance
                    and dists[electrode][chan] <= 2 * gap_tolerance
                ):
                    midpt = self.fill_one_point(
                        final_centroids[electrode][chan_list[i]],
                        final_centroids[electrode][chan_list[i + 1]],
                    )
                    # Shift downstream labels by one
                    midpt_idx = i + 1
                    midpt_id = electrode + str(i + 1)
                    print(f"Interpolating points for channels: " f"{midpt_id}")
                    shift = 1
                    final_centroids[electrode] = self.shift_downstream_labels(
                        electrode, midpt_idx, shift, final_centroids[electrode]
                    )
                    # Insert midpoint into dictionary
                    final_centroids[electrode][midpt_id] = midpt

                # Need to fill in two channels
                elif (
                    dists[electrode][chan] > 2 * gap_tolerance
                    and dists[electrode][chan] <= 3 * gap_tolerance
                ):
                    pt1, pt2 = self.fill_two_points(
                        final_centroids[electrode][chan_list[i]],
                        final_centroids[electrode][chan_list[i + 1]],
                    )
                    # Shift downstream labels by two
                    pt1_idx = i + 1
                    pt1_id = electrode + str(i + 1)
                    pt2_id = electrode + str(i + 2)
                    print(
                        f"Interpolating points for channels: " f"{pt1_id}, " f"{pt2_id}"
                    )
                    shift = 2
                    final_centroids[electrode] = self.shift_downstream_labels(
                        electrode, pt1_idx, shift, final_centroids[electrode]
                    )
                    final_centroids[electrode][pt1_id] = pt1
                    final_centroids[electrode][pt2_id] = pt2
                # Need to fill in three channels
                elif dists[electrode][chan] > 3 * gap_tolerance:
                    pt1, pt2, pt3 = self.fill_three_points(
                        final_centroids[electrode][chan_list[i]],
                        final_centroids[electrode][chan_list[i + 1]],
                    )
                    # Shift downstream labels by three
                    pt1_idx = i + 1
                    pt1_id = electrode + str(i + 1)
                    pt2_id = electrode + str(i + 2)
                    pt3_id = electrode + str(i + 3)
                    print(
                        f"Interpolating points for channels: "
                        f"{pt1_id}, "
                        f"{pt2_id}, "
                        f"{pt3_id}. "
                    )
                    shift = 3
                    final_centroids[electrode] = self.shift_downstream_labels(
                        electrode, pt1_idx, shift, final_centroids[electrode]
                    )
                    final_centroids[electrode][pt1_id] = pt1
                    final_centroids[electrode][pt2_id] = pt2
                    final_centroids[electrode][pt3_id] = pt3
                else:
                    continue
        # Update labeling and ordering for centroids
        final_centroids = self.reassign_labels(final_centroids)
        dists = self.compute_dists(final_centroids)
        return final_centroids, dists

    @classmethod
    def vox_2_xyz(self, final_centroids_voxels, affine):
        """
        Convert finalized dictionary of centroids from CT voxel coordinates
        to xyz coordinates.

        Parameters
        ----------

        final_centroids_voxels: dict()
            Properly labeled dictioanry of centroids in CT voxel coordinates.

        Returns
        -------
            Properly labeled dictionary of centroids in xyz coordinates
        """
        final_centroids_xyz = {}
        for electrode in final_centroids_voxels.keys():
            final_centroids_xyz[electrode] = {}
            for chan in final_centroids_voxels[electrode]:
                final_centroids_xyz[electrode][chan] = apply_affine(
                    affine, final_centroids_voxels[electrode][chan]
                )

        return final_centroids_xyz
