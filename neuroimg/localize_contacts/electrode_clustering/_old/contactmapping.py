# -*- coding: utf-8 -*-

import math
import warnings

import numpy as np
import scipy


class ContactMapping(object):
    """
    Class of functions for mapping contacts.

    Allows mapping to:
    - different regions
    - move electrodes in space
    - find contacts in Euclidean space
    """

    def __init__(self, chanxyz, chanlabels, conn=None):
        self.conn = conn
        self.chanlabels = chanlabels
        self.chanxyz = chanxyz

    def simplest_gain_matrix(self):
        """
        Recompute a new gain matrix based on xyz.

        G = 1 / ( 4*pi * sum(sqrt(( X - X[:, new])^2))^2)

        Returns
        -------
        gain_matrix :
        """
        # NOTE IF YOU MOVE SEEGXYZ ONTO REGXYZ, YOU DIVIDE BY 0, SO THERE IS A PROBLEM
        # reg_xyz = con.centres
        dr = self.conn.centres - self.chanxyz[:, np.newaxis]
        print("Computing simple gain mat!")
        if 0 in dr:
            print(
                "Computing simplest gain matrix will result \
                in error when contact is directly on top of any region!\
                Dividing by 0!"
            )

        ndr = np.sqrt((dr ** 2).sum(axis=-1))
        Vr = 1.0 / (4 * np.pi) / (1 + ndr ** 2)
        return Vr

    def gain_matrix_dipole(
        self, vertices, regmapping, nregions, sensors, orientations, areas
    ):
        """
        Compute gain matrix based on dipole orientations.

        Parameters
        ----------
        vertices             np.ndarray of floats of size n x 3, where n is the number of vertices
        orientations         np.ndarray of floats of size n x 3
        region_mapping       np.ndarray of ints of size n
        sensors              np.ndarray of floats of size m x 3, where m is the number of sensors

        Returns
        -------
        gain_matrix : np.ndarray of size m x n

        """
        nverts = vertices.shape[0]
        nsens = self.chanxyz.shape[0]

        reg_map_mtx = np.zeros((nverts, nregions), dtype=int)
        for i, region in enumerate(regmapping):
            if region >= 0:
                reg_map_mtx[i, region] = 1
        # reg_map_mtx[np.arange(region_mapping.size), region_mapping] = 1.0

        gain_mtx_vert = np.zeros((nsens, nverts))
        for sens_ind in range(nsens):
            a = sensors[sens_ind, :] - vertices
            na = np.sqrt(np.sum(a ** 2, axis=1))
            gain_mtx_vert[sens_ind, :] = (
                areas
                * (np.sum(orientations * a, axis=1) / na ** 3)
                / (4.0 * np.pi * SIGMA)
            )

        return gain_mtx_vert.dot(reg_map_mtx)

    def gain_matrix_inv_square(self, vertices, regmapping, areas):
        """
        Compute a gain matrix using an inverse square fall off (like a mean field model).

        Parameters
        ----------
        vertices             np.ndarray of floats of size n x 3, where n is the number of vertices
        areas                np.ndarray of floats of size n x 3
        region_mapping       np.ndarray of ints of size n
        nregions             int of the number of regions
        sensors              np.ndarray of floats of size m x 3, where m is the number of sensors

        Returns
        -------
        np.ndarray of size m x n
        """
        pass
        nregions = self.conn.region_labels.shape[0]
        nverts = vertices.shape[0]
        nsens = self.chanxyz.shape[0]
        reg_map_mtx = np.zeros((nverts, nregions), dtype=int)
        for i, region in enumerate(regmapping):
            if region >= 0:
                reg_map_mtx[i, region] = 1
        gain_mtx_vert = np.zeros((nsens, nverts))
        for sens_ind in range(nsens):
            a = self.chanxyz[sens_ind, :] - vertices
            na = np.sqrt(np.sum(a ** 2, axis=1))

            # original version
            gain_mtx_vert[sens_ind, :] = areas / (na ** 2)

            # To Do: Refactor to use a more physically accurate way to project source activity
            # adding a 1 in the denominator to softmax the gain matrix
            softmax_inds = np.where(na < 1)[0]
            if len(softmax_inds) > 0:
                # print("na was less than one, so softmaxing here at 1.")
                # epsilon = 1 - a
                # na = np.sqrt(np.sum(a**2, axis=1))
                gain_mtx_vert[sens_ind, softmax_inds] = areas[softmax_inds] / (
                    1 + na[softmax_inds] ** 2
                )

        return gain_mtx_vert.dot(reg_map_mtx)

    def getallcontacts(self, seeg_contact):
        """Get the entire electrode contacts' indices, so that we can modify the corresponding xyz."""
        # get the elec label name
        isleftside = seeg_contact.find("'")
        contacts = []
        for tempcontact in self.chanlabels:
            for idx, s in enumerate(tempcontact):
                if s.isdigit():
                    elec_label = tempcontact[0:idx]
                    break
            contacts.append((elec_label, int(tempcontact[len(elec_label) :])))

        # get indices depending on if it is a left/right hemisphere electrode
        if isleftside != -1:
            elec_label = seeg_contact.split("'")[0]
            electrodeindices = [
                i for i, item in enumerate(self.chanlabels) if elec_label + "'" in item
            ]
        else:
            for idx, s in enumerate(seeg_contact):
                if s.isdigit():
                    elec_label = seeg_contact[0:idx]
                    break
            electrodeindices = [
                i for i, item in enumerate(contacts) if elec_label == item[0]
            ]
        print("\nelec label is %s" % elec_label)
        return electrodeindices

    def _cart2sph(self, x, y, z):
        """
        Transform Cartesian coordinates to spherical.

        Paramters:
        x           (float) X coordinate
        y           (float) Y coordinate
        z           (float) Z coordinate

        :return: radius, elevation, azimuth
        """
        x2_y2 = x ** 2 + y ** 2
        r = math.sqrt(x2_y2 + z ** 2)  # r
        elev = math.atan2(math.sqrt(x2_y2), z)  # Elevation / phi
        az = math.atan2(y, x)  # Azimuth / theta
        return r, elev, az

    def _sph2cart(self, r, elev, az):
        x = r * math.sin(elev) * math.cos(az)
        y = r * math.sin(elev) * math.sin(az)
        z = r * math.cos(elev)
        return x, y, z

    def move_electrode(self, seegind, newloc):
        """
        Move electrode to a new location.

        :param seegind:
        :param newloc:
        :return:
        """
        seeg_contact = self.chanlabels[seegind]
        contact_xyz = self.chanxyz[seegind, :].copy()
        # get all the indices for this electrode
        electrodeindices = self.getallcontacts(seeg_contact=seeg_contact)
        assert len(electrodeindices) > 2

        # get the euclidean distance that will be moved for this electrode
        x_dist = newloc[0] - contact_xyz[0]
        y_dist = newloc[1] - contact_xyz[1]
        z_dist = newloc[2] - contact_xyz[2]
        distancetomove = [x_dist, y_dist, z_dist]
        self.chanxyz[electrodeindices] = self.chanxyz[electrodeindices] + distancetomove

    def findclosestcontact(self, regionind):
        """Find the closest contact to an ezregion."""
        # get the region's xyz coords we want to get
        regionxyz = self.conn.centres[regionind]
        # create a mask of the indices we already moved
        # elec_indices = np.arange(0, self.chanxyz.shape[0])
        # movedmask = [element for i, element in enumerate(elec_indices) \
        #                       if i not in elecmovedindices]
        # create a spatial KD tree -> find closest SEEG contact to region in Euclidean
        # tree = scipy.spatial.KDTree(self.chanxyz[movedmask, :])
        tree = scipy.spatial.KDTree(self.chanxyz)
        near_seeg = tree.query(regionxyz)

        # get the distance and the index at the min
        distance = near_seeg[0]
        seeg_index = near_seeg[1]
        return seeg_index, distance

    def move_electrodetoreg(self, regionind, distance=-1):
        """
        Move the contact and the entire electrode the correct distance.

        Makes it so that the contact is on the ezregion now
        """
        if regionind.size > 1:
            warnings.warn("Need to pass in one region index at a time!")

        if distance == -1:
            warnings.warn("Not moving electrodes, so this call does not do anything!")
            return None, None
        # find the seeg index closest to this region and move it
        seegind, origdistance = self.findclosestcontact(regionind)

        regionxyz = self.conn.centres[regionind, :]
        closest_seegxyz = self.chanxyz[seegind, :].copy()
        seeg_contact = self.chanlabels[seegind]
        # get all the indices for this electrode
        electrodeindices = self.getallcontacts(seeg_contact=seeg_contact)
        print(self.chanlabels[electrodeindices])
        print(self.chanlabels)
        # assert len(electrodeindices) > 2

        # get the euclidean distance that will be moved for this electrode
        x_dist = regionxyz[0] - closest_seegxyz[0]
        y_dist = regionxyz[1] - closest_seegxyz[1]
        z_dist = regionxyz[2] - closest_seegxyz[2]
        distancetomove = [x_dist, y_dist, z_dist]

        # modify the distance in sphereical coordinates
        r, elev, az = self._cart2sph(x_dist, y_dist, z_dist)
        r = r - distance
        x_dist, y_dist, z_dist = self._sph2cart(r, elev, az)
        distancetomove = [x_dist, y_dist, z_dist]

        # createa copy of the seeg_xyz df and modify the electrode
        new_seeg_xyz = self.chanxyz.copy()
        new_seeg_xyz[electrodeindices] = new_seeg_xyz[electrodeindices] + distancetomove

        # modify the object's seeg xyz
        self.chanxyz[electrodeindices] = self.chanxyz[electrodeindices] + distancetomove
        self.gainmat = self.gain_matrix_inv_square()

        print(np.linalg.norm(distancetomove))
        print(origdistance)
        print(distance)
        if distance != -1:
            assert np.linalg.norm(distancetomove) - origdistance <= distance
        return new_seeg_xyz, electrodeindices
