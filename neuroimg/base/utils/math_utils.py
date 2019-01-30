# -*- coding: utf-8 -*-

import math

import numpy as np

from neuroimg.base.config.config import CalculusConfig


class FragilityModel():
    @staticmethod
    def compute_fragilitymetric(minnormpertmat):
        # get dimensions of the pert matrix
        N, T = minnormpertmat.shape
        # assert N < T
        fragilitymat = np.zeros((N, T))
        for icol in range(T):
            fragilitymat[:, icol] = (np.max(minnormpertmat[:, icol]) - minnormpertmat[:, icol]) / \
                                    np.max(minnormpertmat[:, icol])
        return fragilitymat

    @staticmethod
    def compute_minmaxfragilitymetric(minnormpertmat):
        # get dimensions of the pert matrix
        N, T = minnormpertmat.shape
        # assert N < T
        minmax_fragilitymat = np.zeros((N, T))

        # get the min/max for each column in matrix
        minacrosstime = np.min(minnormpertmat, axis=0)
        maxacrosstime = np.max(minnormpertmat, axis=0)

        # normalized data with minmax scaling
        minmax_fragilitymat = -1 * np.true_divide((minnormpertmat - np.matlib.repmat(maxacrosstime, N, 1)),
                                                  np.matlib.repmat(maxacrosstime - minacrosstime, N, 1))
        return minmax_fragilitymat


def normalize_weights(
        weights, percentile=CalculusConfig.WEIGHTS_NORM_PERCENT, remove_diagonal=True, ceil=1.0):
    # Create the normalized connectivity weights:
    if len(weights) > 0:
        normalized_w = np.array(weights)
        if remove_diagonal:
            # Remove diagonal elements
            n_regions = normalized_w.shape[0]
            normalized_w *= 1 - np.eye(n_regions)
        # Normalize with the 95th percentile
        normalized_w = np.array(
            normalized_w /
            np.percentile(
                normalized_w,
                percentile))
        if ceil:
            if ceil is True:
                ceil = 1.0
            normalized_w[normalized_w > ceil] = ceil
        return normalized_w
    else:
        return np.array([])


def compute_in_degree(weights):
    return np.expand_dims(np.sum(weights, axis=1), 1).T


def compute_gain_matrix(locations1, locations2, normalize=95, ceil=1.0):
    n1 = locations1.shape[0]
    n2 = locations2.shape[0]
    projection = np.zeros((n1, n2))
    dist = np.zeros((n1, n2))
    for i1, i2 in product(range(n1), range(n2)):
        dist[i1, i2] = np.abs(
            np.sum((locations1[i1, :] - locations2[i2, :]) ** 2))
        projection[i1, i2] = 1 / dist[i1, i2]
    if normalize:
        projection /= np.percentile(projection, normalize)
    if ceil:
        if ceil is True:
            ceil = 1.0
        projection[projection > ceil] = ceil
    return projection


def cart2sph(x, y, z):
    '''
    Transform Cartesian coordinates to spherical

    Paramters:
    x           (float) X coordinate
    y           (float) Y coordinate
    z           (float) Z coordinate

    :return: radius, elevation, azimuth
    '''
    x2_y2 = x ** 2 + y ** 2
    r = math.sqrt(x2_y2 + z ** 2)  # r
    elev = math.atan2(z, math.sqrt(x2_y2))  # Elevation
    az = math.atan2(y, x)  # Azimuth
    return r, elev, az


def pol2cart(theta, rho):
    '''
    Transform polar coordinates to Cartesian

    Parameters
    theta          (float) angle value
    rho            (float) radius value

    :return: X, Y
    '''
    return rho * math.cos(theta), rho * math.sin(theta)


def azim_proj(pos):
    '''
    Computes the Azimuthal Equidistant Projection of input point in 3D Cartesian Coordinates.
    Imagine a plane being placed against (tangent to) a globe. If
    a light source inside the globe projects the graticule onto
    the plane the result would be a planar, or azimuthal, map
    projection.

    Parameters:
    pos         (list) positions in 3D Cartesian coordinates

    :return: projected coordinates using Azimuthal Equidistant Projection
    '''
    [r, elev, az] = self.cart2sph(pos[0], pos[1], pos[2])
    return self.pol2cart(az, m.pi / 2 - elev)
