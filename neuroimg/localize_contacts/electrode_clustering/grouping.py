import numpy as np
import numpy.linalg as npl
from skimage import measure


class Cluster:
    @classmethod
    def find_clusters(self, maskedCT, pointsOnly=True):
        """
        Function to apply a thresholding based algorithm and then running 
        connected clustering using skimage. This will then return clusters 
        of voxels that belong to distinct contact groups.
        http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label

        Parameters
        ----------
            maskedCT: np.ndarray of 3 dimensions
                Brain-masked CT scan
        
            brainmaskinCT: np.ndarray of 3 dimensions with same shape as masked CT
                Brain mask in CT space
        
            pointsOnly: bool
                True if it is desired to remove background and return 
                a list of all contact coordinates. If False, the function 
                returns a 3D array of the same dimensions as the inputted 
                image array.
        
        Returns
        -------
            allclusters: dict()
                Dictionary of clusters found at various thresholds. The keys are threshold values;
                the values are dictionaries of clusters with their corresponding list of coordinates
                of points in each clsuter.

            numobj: List[int]
                List of the number of clusters found for the range of thresholds
        """
        masktype = "keep"

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
            cluster_labels, numobj[i] = measure.label(
                maskedCT / 255 > currthresh,
                background=0,
                return_num=True,
                connectivity=2,
            )

            print("%.3f \t\t %d" % (currthresh, numobj[i]))

            if pointsOnly is True:
                # filter out all background, i.e. zero-valued, pixels in cluster_labels array
                # tuple of arrays, one for each dim, containing coordinates of nonzero values
                nonzeros = np.nonzero(cluster_labels)
                nonzero_voxs = np.array(
                    list(zip(nonzeros[0], nonzeros[1], nonzeros[2]))
                )

                # go through all identified labels
                for j in range(len(nonzero_voxs)):
                    # Check if coordinate in brainmask is non-zero, i.e. within the brain
                    voxel = nonzero_voxs[j]
                    if (
                        cluster_labels[voxel[0]][voxel[1]][voxel[2]]
                        not in clusters.keys()
                    ):
                        clusters[cluster_labels[voxel[0]][voxel[1]][voxel[2]]] = [
                            nonzero_voxs[j]
                        ]
                    else:
                        clusters[cluster_labels[voxel[0]][voxel[1]][voxel[2]]].append(
                            nonzero_voxs[j]
                        )
                allclusters[float("%.3f" % threshvec[i])] = clusters
            else:
                allclusters[float("%.3f" % threshvec[i])] = cluster_labels
        return allclusters, numobj


class CylindricalGroup:
    @classmethod
    def test_point_in_cylinder(self, pt1, pt2, r, q):
        """
        Tests whether a point q lies within a cylinder with points 
        pt1 and pt2 that define the axis of the cylinder and a 
        specified radius r. Used formulas provided here:
        https://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml

        Parameters
        ----------
        pt1: np.ndarray of dimensions 1 x 3 
            first point to bound the cylinder.
        
        pt2: np.ndarray of dimensions 1 x 3
            second point to bound the cylinder.
        
        r: int
            the radius of the cylinder
        
        q: np.ndarray of dimensions 1 x 3
            the point to test whether it lies within the cylinder
        
        Returns
        -------
            True if q lies in the cylinder, else False.
        """
        vec = pt2 - pt1
        length_sq = npl.norm(vec) ** 2
        radius_sq = r ** 2
        testpt = q - pt1  # set pt1 as origin
        dot = np.dot(testpt, vec)

        # Check if point is within caps of the cylinder
        if dot >= 0 and dot <= length_sq:

            # distance squared to the cylinder axis of the cylinder:
            dist_sq = np.dot(testpt, testpt) - (dot * dot / length_sq)
            if (dist_sq <= radius_sq):
                return True
        return False

    @classmethod
    def points_in_cylinder(self, pt1, pt2, r, points_to_test):
        """
        Tests whether a list of points lies within a cylinder with points pt1 and pt2 that
        define the axis of the cylinder and a specified radius r.

        Parameters
        ----------
            pt1: np.ndarray of dimensions 1 x 3 
                first point to bound the cylinder.
            
            pt2: np.ndarray of dimensions 1 x 3
                second point to bound the cylinder.
            
            r: int
                the radius of the cylinder.

            points_to_test: List[np.ndarray of dimensions 1 x 3]
                list of points to test.
        
        Returns
        -------
            List of all the points that lie in the cylinder
        """
        return [point for point in points_to_test if self.test_point_in_cylinder(pt1, pt2, r, point)]

    @classmethod
    def cylinder_filter(self, elec_in_brain, clusters, radius):
        """
        Applies a cylindrical filtering on the raw threshold-based 
        clustering by generating bounding cylinders given a sparse 
        collection of contact coordinates.
        
        Parameters
        ----------
        elec_in_brain: dict()
            Dictionary of contacts that fall within the brain.
        
        clusters: dict()
            Dictionary of clusters found at a given threshold.
        
        radius: dict()
            Radius with which to form cylindrical boundaries.

        Returns
        -------
            Dictionary of clusters sorted by the cylinder/electrode
            in which they fall. The keys of the dictionary are electrode 
            labels, the values of the dictionary are the cluster points 
            from the threshold-based clustering algorithm that fell 
            into a cylinder.
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
                group_channel_labels[electrode_name] = sorted(
                    group_channel_labels[electrode_name]
                )
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
            sparse_elec_coords[electrode] = [
                elec_in_brain[first],
                elec_in_brain[second],
                elec_in_brain[next_last],
                elec_in_brain[last],
            ]

        # # Hyperparameter: estimating the radius of each electrode in voxels such that
        # # none of the cylinders overlap.
        # radius = radius

        # Perform cylindrical filtering on the points detected by threshold-based clustering
        points_to_test = cluster_points
        clusters_by_cylinder = {}
        for electrode in sparse_elec_coords:
            p1 = sparse_elec_coords[electrode][0]
            p2 = sparse_elec_coords[electrode][-1]
            clusters_by_cylinder[electrode] = {}
            for cluster_id in points_to_test.keys():
                points_list = []
                points_list = self.points_in_cylinder(
                    p1, p2, radius, points_to_test[cluster_id]
                )
                if len(points_list) > 0:
                    clusters_by_cylinder[electrode][cluster_id] = points_list

        return clusters_by_cylinder, sparse_elec_labels, sparse_elec_coords
