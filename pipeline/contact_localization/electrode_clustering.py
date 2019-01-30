import argparse

import nibabel as nb
import numpy as np
import scipy.io
from skimage import measure


def compute_centroids(chanxyzvoxels):
    """
    Function to return the centroids of each channel label given a list of voxels per channel label.

    :param chanxyzvoxels:
    :return:
    """
    pass


def find_electrodes(maskedCT, brainmaskinCT):
    """
    Function to apply a thresholding based algorithm and then running connected clustering using skimage.

    This will then return clusters of voxels that belong to distinct electrode groups.

    http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label

    TODO: test against matlab output and make sure to plot in a jupyter noteobok to ensure that this works.

    :param maskedCT:
    :return:
    """
    masktype = 'keep'

    threshvec = np.arange(0.5, 1, 0.005)

    # create a list to store all objects at each threshold
    numobj = [None] * len(threshvec)
    clusters = []

    # find some optimal threshold
    for i in range(len(threshvec)):
        currthresh = threshvec[i]

        # find connected components in image
        cluster_labels = measure.label(maskedCT / 255 > currthresh,
                                       background=0)

        # go through all identified labels
        for j in range(len(cluster_labels)):
            # get the temporary list of all voxels in cluster
            label = cluster_labels[j]

            # get the thresholded labels that are definitely inside brian mask
            IN = [brainmaskinCT[label] > 0]

            # all outside brain
            if IN == []:
                label = []
            else:
                if masktype == 'keep':
                    # do nothing
                    label = label
                elif masktype == 'partial':
                    # partially crop outside
                    label = label[IN]
                elif masktype == 'remove':
                    # completely remove if outside
                    label = []

            # reset the cluster label if necessary
            cluster_labels[j] = label

        # keep track of clusters per threshold
        clusters.append(cluster_labels)

        # iE = cellfun(@isempty,CC{i}.PixelIdxList)  ;
        # CC{i}.PixelIdxList(iE) = [];
        # CC{i}.NumObjects = length(CC{i}.PixelIdxList) ;
        #
        # NumObj(i) = CC{i}.NumObjects ;


def run_clustering_grouping(maskedCT, threshold, elecinit):
    pass


def main(ctimgfile, brainmaskfile, elecinitfile):
    # load in image volumes using nibabel
    ctimg = nb.load(ctimgfile)
    brainmaskimg = nb.load(brainmaskfile)

    # get brain mask data
    B = brainmaskimg.get_data()

    # extract volume parameters
    nxI = brainmaskimg.header['dims'][0:3]
    dxI = brainmaskimg.header['delta']
    xI = [np.arange(nxi) * dxi - np.mean(np.arange(nxi) * dxi) for nxi, dxi in zip(nxI, dxI)]

    # apply thresholding B
    B[B > 0] = 1

    # mask CT image from skull
    I = ctimg.get_data()
    maskedI = np.multiply(I, B)

    # load inialized electrode
    elecinit = scipy.io.loadmat(elecinitfile)
    elecinitdict = {}
    elecinitpos = elecinit['elecpos']
    elecinitlab = elecinit['label']
    for i in range(len(elecinitpos)):
        elecinitdict[elecinitlab[i]] = elecinitpos[i]

    # run clustering algorithm
    chanxyzvoxels = run_clustering_grouping(maskedI, threshold, elecinitdict)

    # get centroids of each channel
    chanxyz = compute_centroids(chanxyzvoxels)

    return chanxyz, chanxyzvoxels


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ct_nifti_img', help="The CT image volume in its original space.")
    parser.add_argument('brainmask_native_file', help="Brain mask mapped to the CT image space.")
    parser.add_argument('electrode_initialization_file', help="The electrode file with contacts localized to 2 points.")

    parser.add_argument('chanxyz_file', help="The output datafile with all the electrode centroid points labeled.")
    parser.add_argument('clustered_points_file', help="The output datafile with all the electrode points clustered.")
    parser.add_argument('binarized_ct_volume', help='The binarized CT volume.')
    args = parser.parse_args()

    # extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    brainmask_native_file = args.brainmask_native_file
    electrode_initialization_file = args.electrode_initialization_file

    chanxyz_file = args.chanxyz_file
    clustered_points_file = args.clustered_points_file

    threshold = 0.8

    chanxyz, chanxyzvoxels = main(ct_nifti_img, brainmask_native_file, electrode_initialization_file)
