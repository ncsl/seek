import argparse

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from nibabel.affines import apply_affine
from skimage import measure

import sys
import os

sys.path.append("../../")
from neuroimg.base.utils.data_structures_utils import MatReader

sys.path.append("../../../")
from img_pipe.img_pipe import img_pipe

import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import math
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from sklearn.decomposition import PCA


def summary_plots(figurefiledir, final_centroids, elec_in_brain):

    pca_path = os.path.join(figurefiledir, 'summary_PCA.png')
    fig_pca, axs_pca = summary_PCA_plots(pca_path, final_centroids, elec_in_brain)

    l2_path = os.path.join(figurefiledir, 'euclidean_distance_errors.png')
    l2_errors, fig_l2, axs_l2 = l2_error(l2_path, final_centroids, elec_in_brain)

    l2_errorstxt_path = os.path.join(figurefiledir, 'euclidean_distance_errors.txt')
    with open(l2_errorstxt_path, 'w') as f:
        for elec in l2_errors:
            for chan in l2_errors[elec]:
                err = l2_errors[elec][chan]
                f.write("%s %.6f\n" %(chan, err))

    return fig_pca, axs_pca, l2_errors, fig_l2, axs_l2

def summary_PCA_plots(figurefilepath, final_centroids, elec_in_brain):
    if not os.path.exists(os.path.dirname(figurefilepath)):
        os.mkdir(os.path.dirname(figurefilepath))

    sns.set(font_scale=1.1)
    fig, axs = plt.subplots(4, 2, figsize=(10,10))
    axs = axs.flatten()

    grouped_chans_labeled = {}
    for label, contact in elec_in_brain.items():
        if label[1] == "'":
            electrode_name = label[:2] # If electrode name has a '
        else:
            electrode_name = label[0] # If electrode name has no '
        if electrode_name not in grouped_chans_labeled.keys():
            grouped_chans_labeled[electrode_name] = {}
        grouped_chans_labeled[electrode_name][label] = elec_in_brain[label]
    validation = {}
    for electrode in grouped_chans_labeled:
        validation[electrode] = []
        for channel in grouped_chans_labeled[electrode]:
            validation[electrode].append(grouped_chans_labeled[electrode][channel])

    for i, electrode in enumerate(final_centroids):
        # plot stuff
        pca = PCA()
        pca.fit(np.array(list(final_centroids[electrode].values())))
        centroids_pca = pca.transform(np.array(list(final_centroids[electrode].values())))
        centroids_new = pca.inverse_transform(centroids_pca)
        axs[i].scatter(centroids_new[:, 0], np.zeros_like(centroids_new[:, 0]), label='observed')
        axs[i].set_title('PCA Visualization of Electrode %s' % electrode)
        pca.fit(validation[electrode])
        exp_centroids_pca = pca.transform(validation[electrode])
        exp_centroids_new = pca.inverse_transform(exp_centroids_pca)
        axs[i].scatter(exp_centroids_new[:, 0], np.zeros_like(exp_centroids_new[:, 0]), marker='x', c='r', label='expected')
        axs[i].set_title('PCA Validation of Centroids (Electrode: %s)' % electrode)
        axs[i].set_ylim([-0.005, 0.005])
        axs[i].legend()
        axs[i].set_xlabel('PC Coordinates in Voxels along Electrode %s' % electrode)
        for j, chan in enumerate(final_centroids[electrode]):
            axs[i].annotate(chan, (centroids_new[j, 0], 0.0005), size=8.5)

        # set plot paramsx
    fig.tight_layout()

    plt.savefig(figurefilepath,
                box_inches="tight")

    return fig, axs

def l2_error(figurefilepath, final_centroids, elec_in_brain):
    """
    Function that computes the Euclidean distance (l2) between centroids computed by algorithm
    and the validation data, which are centroids manually inputted by user.
    :param final_centroids: dictionary of properly labeled centroids grouped by electrode
    in CT voxels
    :param elec_in_brain: electrode coordinates in CT voxels
    :return error_per_channel: error in each channel (accurate prediction is between 0-3)
    """
    errors_per_channel = {}
    for electrode in final_centroids:
        errors_per_channel[electrode] = {}
        for channel in final_centroids[electrode]:
            observed = final_centroids[electrode][channel] # The coordinates detected by clustering algorithm
            if channel in elec_in_brain:
                expected = elec_in_brain[channel] # Manually labeled validation data
                abs_error = npl.norm(observed-expected)
                errors_per_channel[electrode][channel] = abs_error
            else:
                errors_per_channel[electrode][channel] = float('NaN')


    sns.set(font_scale=1.1)
    fig, axs = plt.subplots(4, 2, figsize=(15,15))
    axs = axs.flatten()

    ymin, ymax = 0, 20
    for i, electrode in enumerate(list(errors_per_channel)):
        y_pos = np.arange(len(errors_per_channel[electrode]))
        axs[i].bar(y_pos, list(errors_per_channel[electrode].values()), align='center', alpha=0.9)
        axs[i].set_xticks(y_pos)
        axs[i].set_xticklabels(list(final_centroids[electrode].keys()))
        axs[i].set_title('Abs. Error By Channel in Electrode %s After Filling Gaps' % electrode)
        axs[i].set_xlabel('Channel')
        axs[i].set_ylabel('Distance')
        axs[i].set_ylim([ymin, ymax])
    fig.tight_layout()

    plt.savefig(figurefilepath,
                box_inches="tight")

    return errors_per_channel, fig, axs



class MaskVolume():
    @classmethod
    def apply_mask(self, mask_data, brain_img_data):
        """
        Applies a binarized brain mask to an input CT image file
        :param mask_data: 3D array of brain mask
        :param brain_img_data: 3D array of CT brain scan
        :return: masked_brain: The masked CT image as a 3D array
        """

        mask_data[mask_data > 0] = 1
        masked_brain = np.multiply(brain_img_data, mask_data)

        return masked_brain

    @classmethod
    def filter_electrodes_bm(self, elec_coords_mm, brainmasked_ct_img):
        """
        Filters out electrodes that do not fall within brain matter of a CT image
        :param elec_coords_mm: A dictionary of contact coordinates in mm space
        :param brainmasked_ct_img: A brainmasked CT NiBabel image object
        :return: elec_in_brain: A dictionary of contact coordinates in CT voxels that fall
            within the binarized brain mask. The keys are individual contact labels, the values
            are the corresponding coordinates in CT space.
        """
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        brainmasked_ct_affine = brainmasked_ct_img.affine
        inv_affine = npl.inv(brainmasked_ct_affine)

        # Convert contact xyz coordinates to CT voxels
        elec_coords_CTvox = {}
        for label, contact in elec_coords_mm.items():
            elec_coords_CTvox[label] = apply_affine(inv_affine, contact)

        # Filter out electrodes not within brain mask
        elec_in_brain = {}
        for label, contact in elec_coords_CTvox.items():
            if brainmasked_ct_data[int(contact[0]), int(contact[1]), int(contact[2])] != 0:
                elec_in_brain[label] = contact
        return elec_in_brain

    @classmethod
    def sort_contacts(self, elec_in_brain):
        """
        Groups the individual contacts by the electrode to which they correspond
        :param elec_in_brain: A dictionary of contact coordinates in CT voxels that fall
            within the brain matter.
        :return: voxels_per_electrode: A dictionary of contact coordinates in CT voxels that
            fall within the brain. The keys are electrode labels, and the values are lists of
            the coordinates of all the contacts that correspond to a given electrode.
        """
        voxels_per_electrode = {}
        for label, contact in elec_in_brain.items():
            if label[1] == "'":
                electrode_name = label[:2]  # If electrode name has a '
            else:
                electrode_name = label[0]  # If electrode name has no '
            if electrode_name in voxels_per_electrode.keys():
                voxels_per_electrode[electrode_name].append(contact)
            else:
                voxels_per_electrode[electrode_name] = [contact]


class Cluster():
    @classmethod
    def find_clusters(self, maskedCT, pointsOnly=True):
        """
        Function to apply a thresholding based algorithm and then running connected clustering using skimage.
        This will then return clusters of voxels that belong to distinct contact groups.
        http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label
        :param maskedCT: Brain-masked CT scan
        :param brainmaskinCT: Brain mask in CT space
        :param pointsOnly: True if it is desired to remove background and return a list of all contact coordinates.
            If False, the function returns a 3D array of the same dimensions as the inputted image array.
        :return:
            allclusters: Dictionary of clusters found at various thresholds. The keys are threshold values;
                the values are dictionaries of clusters with their corresponding list of coordinates
                of points in each clsuter.
            numobj: List of the number of clusters found for the range of thresholds
        """
        masktype = 'keep'

        threshvec = np.arange(0.63, 0.64, 0.005)
        # create a list to store number of clusters at each threshold
        numobj = [None] * len(threshvec)
        allclusters = {}

        print("Running clustering with various thresholds:")
        print("Threshold \t Num Clusters Found")

        # find some optimal threshold
        for i in range(len(threshvec)):
            clusters = {}
            currthresh = threshvec[i]

            # find connected components in image
            cluster_labels, numobj[i] = measure.label(maskedCT / 255 > currthresh,
                                                      background=0, return_num=True, connectivity=2)

            print("%.3f \t\t %d" % (currthresh, numobj[i]))

            if pointsOnly is True:
                # filter out all background, i.e. zero-valued, pixels in cluster_labels array
                # tuple of arrays, one for each dim, containing coordinates of nonzero values
                nonzeros = np.nonzero(cluster_labels)
                nonzero_voxs = np.array(list(zip(nonzeros[0], nonzeros[1], nonzeros[2])))

                # go through all identified labels
                for j in range(len(nonzero_voxs)):
                    # Check if coordinate in brainmask is non-zero, i.e. within the brain
                    voxel = nonzero_voxs[j]
                    if (cluster_labels[voxel[0]][voxel[1]][voxel[2]] not in clusters.keys()):
                        clusters[cluster_labels[voxel[0]][voxel[1]][voxel[2]]] = [nonzero_voxs[j]]
                    else:
                        clusters[cluster_labels[voxel[0]][voxel[1]][voxel[2]]].append(nonzero_voxs[j])
                allclusters[float("%.3f" % threshvec[i])] = clusters
            else:
                allclusters[float("%.3f" % threshvec[i])] = cluster_labels
        return allclusters, numobj


class CylindricalGroup():
    @classmethod
    def test_point_in_cylinder(self, pt1, pt2, r, q):
        """
        Tests whether a point q lies within a cylinder with points pt1 and pt2 that
        define the axis of the cylinder and a specified radius r. Used formulas provided here:
        https://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml
        :param pt1: first point to bound the cylinder
        :param pt2: second point to bound the cylinder
        :param r: the radius of the cylinder
        :param q: the point to test whether it lies within the cylinder
        :return True if q lies in the cylinder, else False.
        """
        vec = pt2 - pt1
        length_sq = npl.norm(vec) ** 2
        radius_sq = r ** 2
        testpt = q - pt1  # set pt1 as origin
        dot = np.dot(testpt, vec)

        # Check if point is within caps of the cylinder
        if (dot >= 0 and dot <= length_sq):

            # distance squared to the cylinder axis of the cylinder:
            dist_sq = np.dot(testpt, testpt) - (dot * dot / length_sq)
            if (dist_sq <= r):
                return True
        return False

    @classmethod
    def points_in_cylinder(self, pt1, pt2, r, points_to_test):
        """
        Tests whether a list of points lies within a cylinder with points pt1 and pt2 that
        define the axis of the cylinder and a specified radius r.
        :param pt1: first point to bound the cylinder
        :param pt2: second point to bound the cylinder
        :param r: the radius of the cylinder
        :param points_to_test: list of points to test
        :return: points_list: all the points that lie in the cylinder
        """
        points_list = []
        for point in points_to_test:
            if self.test_point_in_cylinder(pt1, pt2, r, point) == True:
                points_list.append(point)
        return points_list

    @classmethod
    def cylinder_filter(self, elec_in_brain, clusters, r):
        """
        Applies a cylindrical filtering on the raw threshold-based clustering by generating
        bounding cylinders given a sparse collection of contact coordinates.
        :param elec_in_brain: Dictionary of contacts that fall within the brain
        :param clusters: Dictionary of clusters found at a given threshold
        :param radius: Radius with which to form cylindrical boundaries
        :return: clusters_by_cylinder: Dictionary of clusters sorted by the cylinder/electrode
            in which they fall. The keys of the dictionary are electrode labels, the values of the
            dictionary are the cluster points from the threshold-based clustering algorithm that
            fell into a cylinder.
        """

        # Separate components of all points found in clustering algorithm
        # and group them.
        cluster_points = {}
        cluster_idata = []
        cluster_jdata = []
        cluster_kdata = []

        for cluster_id in clusters:
            cluster_points[cluster_id] = []
            for point in clusters[cluster_id]:
                cluster_points[cluster_id].append(point)
                cluster_idata.append(point[0])
                cluster_jdata.append(point[1])
                cluster_kdata.append(point[2])

        # Group all the contacts by the electrode they correspond to.
        group_channel_labels = {}
        for label in elec_in_brain:
            if label[1] == "'":
                electrode_name = label[:2]  # If electrode name has a '
                chan_num = label[2:]
            else:
                electrode_name = label[0]  # If electrode name has no '
                chan_num = label[1:]
            if electrode_name in group_channel_labels.keys():
                group_channel_labels[electrode_name].append(int(chan_num))
                group_channel_labels[electrode_name] = sorted(group_channel_labels[electrode_name])
            else:
                group_channel_labels[electrode_name] = [int(chan_num)]

        sparse_elec_labels = {}
        sparse_elec_coords = {}

        # In case the channel labels were not in order, we now sort them in order
        # and then consider only the most medial and surface-level contacts.
        for electrode in group_channel_labels:
            labels = group_channel_labels[electrode]
            first_num = labels[0]
            second_num = labels[1]
            last_num = labels[-1]
            next_last_num = labels[-2]
            first = "%s%d" % (electrode, first_num)
            second = "%s%d" % (electrode, second_num)
            next_last = "%s%d" % (electrode, next_last_num)
            last = "%s%d" % (electrode, last_num)
            sparse_elec_labels[electrode] = [first, second, next_last, last]
            sparse_elec_coords[electrode] = [elec_in_brain[first], elec_in_brain[second],
                                             elec_in_brain[next_last], elec_in_brain[last]]

        # Hyperparameter: estimating the radius of each electrode in voxels such that
        # none of the cylinders overlap.
        radius = r

        # Perform cylindrical filtering on the points detected by threshold-based clustering
        points_to_test = cluster_points
        clusters_by_cylinder = {}
        for electrode in sparse_elec_coords:
            p1 = sparse_elec_coords[electrode][0]
            p2 = sparse_elec_coords[electrode][-1]
            clusters_by_cylinder[electrode] = {}
            for cluster_id in points_to_test.keys():
                points_list = []
                points_list = self.points_in_cylinder(p1, p2, radius, points_to_test[cluster_id])
                if len(points_list) > 0:
                    clusters_by_cylinder[electrode][cluster_id] = points_list

        return clusters_by_cylinder, sparse_elec_labels, sparse_elec_coords


class PostProcessor():
    @classmethod
    def process_abnormal_clusters(self, clusters, elec_in_brain, cylindrical_filter_clusters, sparse_elec_labels):
        """
        Truncates the clusters closest to the skull, which tend to be oversized, and
        separates clusters that appear to have grouped two contacts together. The cluster is filtered
        with 90% quantile filtering, using a coordinate from user-input as the mean for the cluster.
        :param clusters: Dictionary of clusters found at a given threshold
        :param elec_in_brain: Dictionary of all the contacts found within the brain
        :param cylindrical_filter_clusters: Dictionary of points first sorted by the cylindrical electrode
            in which they fall and then the cluster in which they fall. The keys of the dictionary are
            electrode labels, the values of the dictionary are the cluster points from the threshold-based
            clustering algorithm that fell into a cylinder.
        :param sparse_elec_labels: sparse dictionary of electrode labels given by user input
        :return: clusters: Updated dictionary of clusters found along each electrode
            with the skull clusters having been resized through quantile filtering and large clusters having been separated.

        """
        # Group cylindrically bounded clusters by electrodes
        clusters_in_cylinder_by_elec = {}
        for electrode in cylindrical_filter_clusters:
            clusters_in_cylinder_by_elec[electrode] = []
            for cluster_id in cylindrical_filter_clusters[electrode]:
                for i in range(len(cylindrical_filter_clusters[electrode][cluster_id])):
                    point = cylindrical_filter_clusters[electrode][cluster_id][i]
                    clusters_in_cylinder_by_elec[electrode].append(point)

        cluster_labels_per_electrode = {}
        for electrode in clusters_in_cylinder_by_elec:
            cluster_labels_per_electrode[electrode] = []
            for cluster_id in clusters:
                a = map(tuple, clusters[cluster_id]) # Need to do this to convert to set and intersect in O(1)
                b = map(tuple, clusters_in_cylinder_by_elec[electrode])
                intxn = list(set(a) & set(b))
                if len(intxn) > 0:
                    cluster_labels_per_electrode[electrode].append(cluster_id)

        cluster_sizes_by_electrode = {}
        for electrode in cluster_labels_per_electrode:
            cluster_sizes_by_electrode[electrode] = []
            for cluster_id in cluster_labels_per_electrode[electrode]:
                cluster = clusters[cluster_id]
                cluster_sizes_by_electrode[electrode].append(len(cluster))

        # Identify abnormal clusters
        merged_clusters = []
        skull_clusters = []
        for electrode, cluster_sizes in cluster_sizes_by_electrode.items():
            for i in range(len(cluster_sizes)):
                # Hyperparameters: threshold for which we consider a cluster to be a skull cluster or fused cluster
                if cluster_sizes[i] > 200:
                    skull_clusters.append(cluster_labels_per_electrode[electrode][i])
                elif cluster_sizes[i] >= 50 and cluster_sizes[i] <= 200:
                    merged_clusters.append(cluster_labels_per_electrode[electrode][i])

        result = {}
        for electrode in cylindrical_filter_clusters:
            result[electrode] = {}
            for cluster_id in cylindrical_filter_clusters[electrode]:
                result[electrode][cluster_id] = []
                for point in cylindrical_filter_clusters[electrode][cluster_id]:
                    result[electrode][cluster_id].append(point)

        resized_clusters = {}
        # Hyperparameter: qtile
        qtile = 0.875
        for cluster_id in skull_clusters:
            cluster_points = clusters[cluster_id]
            var = np.var(cluster_points, axis=0)
            for electrode in cluster_labels_per_electrode:
                if cluster_id in cluster_labels_per_electrode[electrode]:
                    last_contact_label = sparse_elec_labels[electrode][-1]
                    last_contact = elec_in_brain[last_contact_label]
                    resized_clusters[cluster_id] = []
                    for point in cluster_points:
                        z = npl.norm(np.divide((point - last_contact), np.sqrt(var)))
                        # Only accept up to 90% quantile
                        if norm.cdf(z) <= qtile:
                            resized_clusters[cluster_id].append(point)

        # Update resulting dictionary after resizing skull clusters
        for electrode in result:
            for cluster_id in resized_clusters:
                if cluster_id in result[electrode].keys():
                    result[electrode][cluster_id] = resized_clusters[cluster_id]

        # Separate fused clusters into two new clusters and store them in a dictionary
        unmerged_clusters = {}
        max_cluster_id = max(clusters.keys())
        for cluster_id in merged_clusters:
            pts = clusters[cluster_id]
            var = np.var(pts, axis=0)
            var = npl.norm(var)
            for electrode in cluster_labels_per_electrode:
                if cluster_id in cluster_labels_per_electrode[electrode]:
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
                            unmerged_clusters[electrode][max_cluster_id + 1].append(point)
                        if norm.cdf(z2) <= qtile:
                            unmerged_clusters[electrode][max_cluster_id + 2].append(point)
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
                    result[electrode][cluster_id] = unmerged_clusters[electrode][cluster_id]

        relabeled_clusters = {}
        for electrode in result:
            num=1
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
        :param: p1: first point to compute midpoint.
        :param: p2: second point to compute midpoint.
        :return: midpoint of the two specified points.
        """
        return (p1+p2)/2

    @classmethod
    def fill_two_points(self, p1, p2):
        """
        Computes trisection points.
        :param: p1: first point to compute trisection points.
        :param: p2: second point to compute trisection points.
        :return: two points which trisect the line segment formed by p1 and p2
        """
        return (2/3)*p1+(1/3)*p2, (1/3)*p1+(2/3)*p2

    @classmethod
    def fill_three_points(self, p1, p2):
        """
        Computes quadrisection points.
        :param: p1: first point to compute quadrisection points.
        :param: p2: second point to compute qaudrisection points.
        :return: three points which quadrisect the line segment formed by p1 and p2
        """
        return (1/4)*p1 + (3/4)*p2, self.fill_one_point(p1,p2), (3/4)*p1 + (1/4)*p2

    @classmethod
    def shift_downstream_labels(self, electrode_name, idx, shift_factor, chan_dict):
        """
        Helps to relabel downstream labels from a given index to maintain sorted order.
        :param: electrode_name: string that contains electrode name (e.g. L')
        :param: idx: index from which shifting should start
        :param: shift_factor: how much to shift downstream labels by
        :param: chan_dict: dictionary of channels with standard labeling and their
        corresponding coordinates.
        :return: updated version of chan_dict with shifted labels to allow for easy insertion
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
        Compute the distances between a given channel and its immediate neighbor on the right.
        The last channel has a distance set to 0.
        :param: centroid_dict: a dictionary with keys being electrode names and values being
        dictionaries consisting of entries of channel names and their corresponding centroid
        coordinates.
        :return: dists: a dictionary with keys being electrode names and values being dictionaries
        consisting of entries of channel names and their corresponding distances as described
        above.
        """

        dists = {elec: {} for elec in list(centroid_dict.keys())}
        for electrode in centroid_dict:
            channels = list(centroid_dict[electrode].keys())
            num_chans = len(channels)
            centrs = centroid_dict[electrode]
            for i in range(num_chans - 1):
                cur = np.array(centrs[channels[i]])
                nxt = np.array(centrs[channels[i+1]])
                dists[electrode][channels[i]] = npl.norm(cur - nxt)
            dists[electrode][channels[-1]] = 0
        return dists

    @classmethod
    def reassign_labels(self, centroid_dict):
        """
        Assigns labels which follow the standard labeling convention for SEEG electrodes.
        :param centroid_dict: A dictionary where the keys are each electrode name and values
        are dictionaries with entries of channels (does not necessarily have to follow standard
        labeling convention) and their corresponding coordinates.
        :return: result: A dictionary where channels are correctly assigned and sorted in order.
        If the electrode name has a ', then the first entry in the dictionary will be the smallest
        number (usually 1) and every subsequent entry will be larger. Otherwise, the first entry
        in the dictionary will be the largest number and every subsequent entry will be smaller.
        """
        pca = PCA()
        result = {}
        for electrode in centroid_dict:
            result[electrode] = {}
            pca.fit(np.array(list(centroid_dict[electrode].values())))
            centroids_pca = pca.transform(np.array(list(centroid_dict[electrode].values())))
            centroids_new = pca.inverse_transform(centroids_pca)
            sorted_idxs = np.argsort(centroids_new[:, 0])
            sorted_orig_ids = np.array(list(centroid_dict[electrode].keys()))[sorted_idxs]
            # Left side of the brain
            assert len(sorted_idxs) == len(centroid_dict[electrode].keys())
            if electrode[-1] == "'":
                for i, chan in enumerate(sorted_idxs):
                    new_id = electrode + str(i+1)
                    result[electrode][new_id] = centroid_dict[electrode][sorted_orig_ids[i]]
            # Right side of the brain
            else:
                for i, chan in enumerate(sorted_idxs):
                    new_id = electrode + str(len(sorted_orig_ids)-i)
                    result[electrode][new_id] = centroid_dict[electrode][sorted_orig_ids[i]]
        return result

    @classmethod
    def fill_gaps(self, final_centroids, dists, gap_tolerance):
        """
        Assist in filling in gaps in clustering.
        :param: final_centroids: a dictionary with keys being electrode names and values being
        dictionaries consisting of entries of channel names and their corresponding centroid
        coordinates.
        :param: dists: a dictionary with keys being electrode names and values being dictionaries
        consisting of entries of channel names and their corresponding distance to their most
        immediate neighbor on their right. The rightmost channel for a given electrode has a
        distance set to 0.
        :param: gap_tolerance: maximum distance we will allow two adjacent contacts to be before
        no longer considering them adjacent.
        :return: updated versions of final_centroids with adjusted labeling and dists with updated
        distances.
        """
        for electrode in dists:
            for i, chan in enumerate(dists[electrode]):
                if i == len(list(dists[electrode].keys())) - 1:
                    continue
                chan_list = np.array(list(dists[electrode].keys()))
                # Hyperparameters: What distance to consider having skipped
                # one, two, or three channels
                # Need to fill in one channel
                if dists[electrode][chan] > gap_tolerance and dists[electrode][chan] <= 2*gap_tolerance:
                    midpt = self.fill_one_point(final_centroids[electrode][chan_list[i]], final_centroids[electrode][chan_list[i+1]])
                    midpt_idx = i+1
                    midpt_id = electrode + str(i+1)
                    # Shift downstream labels by one
                    shift = 1
                    final_centroids[electrode] = self.shift_downstream_labels(electrode, midpt_idx, shift, final_centroids[electrode])
                    # Insert midpoint into dictionary
                    final_centroids[electrode][midpt_id] = midpt

                # Need to fill in two channels
                elif dists[electrode][chan] > 2*gap_tolerance and dists[electrode][chan] <= 3*gap_tolerance:
                    pt1, pt2 = self.fill_two_points(final_centroids[electrode][chan_list[i]], final_centroids[electrode][chan_list[i+1]])
                    pt1_idx = i+1
                    pt2_idx = i+2
                    pt1_id = electrode + str(i+1)
                    pt2_id = electrode + str(i+2)
                    # Shift downstream labels by two
                    shift = 2
                    final_centroids[electrode] = self.shift_downstream_labels(electrode, pt1_idx, shift, final_centroids[electrode])
                    final_centroids[electrode][pt1_id] = pt1
                    final_centroids[electrode][pt2_id] = pt2
                # Need to fill in three channels
                elif dists[electrode][chan] > 3*gap_tolerance:
                    pt1, pt2, pt3 = self.fill_three_points(final_centroids[electrode][chan_list[i]], final_centroids[electrode][chan_list[i+1]])
                    # Shift downstream labels by three
                    pt1_idx = i+1
                    pt2_idx = i+2
                    pt3_idx = i+3
                    pt1_id = electrode + str(i+1)
                    pt2_id = electrode + str(i+2)
                    pt3_id = electrode + str(i+3)
                    shift = 3
                    final_centroids[electrode] = self.shift_downstream_labels(electrode, pt1_idx, shift, final_centroids[electrode])
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
    def vox_2_xyz(self, final_centroids_voxels):
        """
        Convert finalized dictionary of centroids from CT voxel coordinates
        to xyz coordinates.
        :param final_centroids_voxels: Properly labeled dictioanry of centroids
        in CT voxel coordinates
        :return final_centroids_xyz: Properly labeled dictionary of centroids in
        xyz coordinates
        """
        final_centroids_xyz = {}
        for electrode in final_centroids_voxels.keys():
            final_centroids_xyz[electrode] = {}
            for chan in final_centroids_voxels[electrode]:
                final_centroids_xyz[electrode][chan] = final_centroids_voxels[electrode][chan]
        return final_centroids_xyz



def load_data(ct_scan, brainmask_ct, elecinitfile):
    """
    Load each brain image scan as a NiBabel image object
    :param ct_scan: Nifti image file of CT scan
    :param brainmask_ct: Nifti image file of corresponding brain mask in CT voxels
    :param elecinitfile: Space-delimited text file of contact labels and contact
        coordinates in mm space
    :return:
        ct_img: NiBabel image object of CT scan input
        brainmask_ct: NiBabel image object of brain mask in CT
        elecinitfile: A dictionary of contact coordinates in mm space. Keys are
            individual contact labels, and values are the corresponding coordinates
            in mm space.
    """

    ct_img = nb.load(ct_scan)
    bm_ct_img = nb.load(brainmask_ct)
    elec_coords_mm = {}
    with open(elecinitfile) as f:
        for l in f:
            row = l.split()
            elec_coords_mm[row[0]] = np.array([float(row[1]), float(row[2]), float(row[3])])
    return ct_img, bm_ct_img, elec_coords_mm


def compute_centroids(chanxyzvoxels):
    """
    Function to return the centroids of each channel label given a list of voxels per channel label.

    :param chanxyzvoxels: dictionary of electrodes and corresponding dictionary
    of channels to centroid xyz coordinates.
    :return: centroids: list of centroids
    """
    centroids = {}
    for channel, voxels in chanxyzvoxels.items():
        i = j = k = count = 0
        for vox in voxels:
            i += vox[0]
            j += vox[1]
            k += vox[2]
            count += 1
        centroids[channel] = np.array([int(i / count), int(j / count), int(k / count)])
    return centroids


def atlas_and_mask(final_centroids_xyz, PATIENT_OUTPUT_DIR):
    """
    Map centroids to an atlas (e.g. Desikan-Killiany, Destriuex) and apply
    white matter and brain masks to label centroids as white matter or out of
    the brain.
    :param final_centroids_xyz: dictionary of electrodes and corresponding
    dictionary of channels to centroid xyz coordinates.
    :return elec_labels_destriuex: array of contacts labeled with Destriuex atlas
    :return elec_labels_DKT: array of contacts labeled with Desikan-Killiany atlas
    """

    patid = os.path.basename(os.path.normpath(PATIENT_OUTPUT_DIR))

    elecs_dir = os.path.join(PATIENT_OUTPUT_DIR, 'elecs')
    if not os.path.exists(elecs_dir):
        os.mkdir(elecs_dir)

    fsoutmatnamexyz = os.path.join(elecs_dir, '%s_clusteredelec_xyz.mat' % (patid))

    # Save centroids as .mat file with attributes eleclabels, which stores
    # channel name, electrode type, and depth/grid/strip/other and elecmatrix,
    # which stores the centroid coordinates
    eleclabels = []
    elecmatrix = []
    for elec in final_centroids_xyz:
        for chan in final_centroids_xyz[elec]:
            label = [[chan], ['stereo'], ['depth']]
            eleclabels.append(label)
            elecmatrix.append(final_centroids_xyz[elec][chan])
    mat = {'eleclabels': eleclabels, 'elecmatrix': elecmatrix}
    scipy.io.savemat(fsoutmatnamexyz, mat)


    # Apply Atlases, white matter mask, and brainmask
    freeCoG = img_pipe.freeCoG(subj=patid, hem="lh")
    freeCoG.convert_fsmesh2mlab()
    elec_labels_destriuex = freeCoG.label_elecs(elecfile_prefix="la04_clusteredelec_xyz",
                                    atlas_depth="destriuex")
    temp = scipy.io.loadmat(fsoutmatnamexyz)
    destriuexname = os.path.join(elecs_dir, '%s_clustered_elec_xyz_destriuex.mat' % (patid))
    scipy.io.savemat(destriuexname, temp)

    elec_labels_DKT = freeCoG.label_elecs(elecfile_prefix="la04_clusteredelec_xyz",
                                    atlas_depth="desikan-killiany")
    temp = scipy.io.loadmat(fsoutmatnamexyz)
    dktname = os.path.join(elecs_dir, '%s_clustered_elec_xyz_DKT' % (patid))
    scipy.io.savemat(dktname, temp)

    elec_labels_destriuex = wm_and_brainmask(final_centroids_xyz, destriuexname, PATIENT_OUTPUT_DIR)
    elec_labels_DKT = wm_and_brainmask(final_centroids_xyz, dktname, PATIENT_OUTPUT_DIR)

    return elec_labels_destriuex, elec_labels_DKT

def wm_and_brainmask(final_centroids_xyz, elecfile, PATIENT_OUTPUT_DIR):
    """
    Apply white matter and brainmask labels to final centroid output and save
    in .mat files
    :param final_centroids_xyz:
    :param elecfile:
    :param PATIENT_OUTPUT_DIR:
    :return elec_labels:
    """
    dat = scipy.io.loadmat(elecfile)
    elecmatrix = dat['elecmatrix']
    anatomy_orig = dat['anatomy']
    eleclabels = dat['eleclabels']

    # Load white matter and brain masks
    wmpath = os.path.join(PATIENT_OUTPUT_DIR, 'CT', 'wm_ct.nii.gz')
    bmpath = os.path.join(PATIENT_OUTPUT_DIR, 'CT', 'brainmask_native_inct.nii.gz')
    wm_img = nb.load(wmpath)
    wm_dat = wm_img.get_data()
    bm_img = nb.load(bmpath)
    bm_dat = bm_img.get_data()

    affine = npl.inv(bm_img.affine)

    wm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    bm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    anatomy = np.zeros((anatomy_orig.shape[0], anatomy_orig.shape[1] + 2), dtype=object)
    for i, label in enumerate(anatomy_orig):
        chan = str(label[0][0]).strip()
        for elec in final_centroids_xyz:
            if chan in final_centroids_xyz[elec].keys():
                pt = apply_affine(affine, final_centroids_xyz[elec][chan])
                wm_label[i] = wm_dat[int(pt[0]), int(pt[1]), int(pt[2])] > 0
                bm_label[i] = bm_dat[int(pt[0]), int(pt[1]), int(pt[2])] > 0
    anatomy[:, 0:anatomy_orig.shape[1]] = anatomy_orig
    anatomy[:, anatomy_orig.shape[1]] = wm_label
    anatomy[:, anatomy_orig.shape[1] + 1] = bm_label

    save_dict = {'elecmatrix': elecmatrix, 'anatomy': anatomy, 'eleclabels': eleclabels}
    scipy.io.savemat(elecfile, save_dict)

    elec_labels = anatomy

    return elec_labels


# <<<<<<< HEAD
#     # print(brainmaskimg.header)
#
#     # extract volume parameters
#     dimattribname = "dim"
#     pixdim_attribname = 'pixdim'
#
#     nxI = brainmaskimg.header[dimattribname][0:3]
#     dxI = brainmaskimg.header[pixdim_attribname]
#     xI = [np.arange(nxi) * dxi - np.mean(np.arange(nxi) * dxi) for nxi, dxi in zip(nxI, dxI)]
#
#     # load inialized electrode
#     matreader = MatReader()
#     elecinit = matreader.loadmat(elecinitfile)
#     elecinit = elecinit['elecf']
#     # print(elecinit)
#     print(elecinit.keys())
#
#
#     elecinitdict = {}
#     elecinitpos = elecinit['elecpos']
#     elecinitlab = elecinit['label']
#     for i in range(len(elecinitpos)):
#         elecinitdict[elecinitlab[i]] = elecinitpos[i]
#
def main(ctimgfile, brainmaskfile, elecinitfile):
    # Load data
    ct_img, bm_ct_img, elec_coords_mm = load_data(ctimgfile, brainmaskfile, elecinitfile)
    ct_data = ct_img.get_fdata()
    bm_ct_data = bm_ct_img.get_fdata()

    # define pipeline objects to run algorithm
    maskpipe = MaskVolume()
    clusterpipe = Cluster()
    grouppipe = CylindricalGroup()
    postprocesspipe = PostProcessor()

    # Apply masking
    brainmasked_ct_data = maskpipe.apply_mask(bm_ct_data, ct_data)
    brainmasked_ct_img = nb.Nifti1Image(brainmasked_ct_data, bm_ct_img.affine)

    # Filtering out electrodes not within brainmask
    elec_in_brain = maskpipe.filter_electrodes_bm(elec_coords_mm, brainmasked_ct_img)
    voxels_per_electrode = maskpipe.sort_contacts(elec_in_brain)

    # Runs threshold-based clustering algorithm over brainmasked CT
    # for thresholds between 0.62 and 0.65 with a step of 0.005
    clusters, numobj = np.array(clusterpipe.find_clusters(brainmasked_ct_data))

    # Threshold for applying naive clustering algorithm
    threshold = 0.630

    # Cluster by cylinder
    radius = 15
    clusters_by_cylinder, sparse_elec_labels, sparse_elec_coords = grouppipe.cylinder_filter(elec_in_brain, clusters[threshold], radius)

    # Begin postprocessing steps
    processed_clusters = postprocesspipe.process_abnormal_clusters(clusters[threshold], elec_in_brain, clusters_by_cylinder, sparse_elec_labels)

    # Compute centroids for filling gaps
    gap_tolerance = 12.5 # maximum distance between two adjacent nodes
    centroids = {}
    for electrode in processed_clusters:
        centroids[electrode] = compute_centroids(processed_clusters[electrode])
    centroids = postprocesspipe.reassign_labels(centroids)

    dists = postprocesspipe.compute_dists(centroids)
    final_centroids_voxels, dists = postprocesspipe.fill_gaps(centroids, dists, gap_tolerance)
    final_centroids_xyz = postprocesspipe.vox_2_xyz(final_centroids_voxels)

    summary_plots_path = "/Users/ChesterHuynh/sarmalab/contact_localization_proj/fs_outputdata/figures/"
    fig_pca, axs_pca, l2_errors, fig_l2, axs_l2 = summary_plots(summary_plots_path, final_centroids_voxels, elec_in_brain)

    # Sum over all individual errors
    total_err = 0
    for elec in l2_errors:
        for chan in l2_errors[elec]:
            if not math.isnan(l2_errors[elec][chan]):
                total_err += l2_errors[elec][chan]

    # Output labeled .mat files with atlas, white matter, and brainmask information
    PATIENT_OUTPUT_DIR = "/Users/ChesterHuynh/sarmalab/contact_localization_proj/fs_outputdata/outputfiles/la04"
    atlas_and_mask(final_centroids_xyz, PATIENT_OUTPUT_DIR)

    return final_centroids_voxels, final_centroids_xyz, brainmasked_ct_img


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ct_nifti_img', help="The CT image volume in its original space.")
    parser.add_argument('brainmask_native_file', help="Brain mask mapped to the CT image space.")
    parser.add_argument('electrode_initialization_file', help="The electrode file with contacts localized to 2 points.")
    parser.add_argument('chanxyz_file', help="The output datafile with all the electrode centroid points labeled.")
    parser.add_argument('clustered_points_file', help="The output datafile with all the electrode points clustered.")
    parser.add_argument("clustered_voxels_file", help="the voxels output datafile")
    parser.add_argument('binarized_ct_volume', help='The binarized CT volume.')
    args = parser.parse_args()

    # extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    brainmask_native_file = args.brainmask_native_file
    electrode_initialization_file = args.electrode_initialization_file
    chanxyz_file = args.chanxyz_file
    clustered_points_file = args.clustered_points_file
    clustered_voxels_file = args.clustered_voxels_file
    binarized_ct_file = args.binarized_ct_volume

    final_centroids_voxels, final_centroids_xyz, binarized_ct_img = main(ct_nifti_img, brainmask_native_file, electrode_initialization_file)

    #save output files
    with open(clustered_points_file, 'w') as f:
        for elec in final_centroids_xyz:
            for chan in final_centroids_xyz[elec]:
                pt = final_centroids_xyz[elec][chan]
                f.write("%s %.6f %.6f %.6f\n" % (chan, pt[0], pt[1], pt[2]))

    with open(clustered_voxels_file, 'w') as f:
        for elec in final_centroids_voxels:
            for chan in final_centroids_voxels[elec]:
                vox = final_centroids_voxels[elec][chan]
                f.write("%s %.6f %.6f %.6f\n" % (chan, vox[0], vox[1], vox[2]))

    binarized_ct_img.to_filename(binarized_ct_file)
