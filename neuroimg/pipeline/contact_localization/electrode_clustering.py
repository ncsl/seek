import argparse
import os
import sys

import nibabel as nb
import numpy as np
import numpy.linalg as npl
import scipy.io
from nibabel.affines import apply_affine

sys.path.append("../../../")

from neuroimg.base.utils import MatReader
from neuroimg.localize_contacts.electrode_clustering.mask import MaskVolume
from neuroimg.localize_contacts.electrode_clustering.grouping import Cluster, CylindricalGroup
from neuroimg.localize_contacts.electrode_clustering.postprocess import PostProcessor
from neuroimg.localize_contacts.freecog_labeling.utils import convert_fsmesh2mlab, label_elecs


def load_data(ct_scan, brainmask_ct=None):
    """
    Load brain image scans as NiBabel image objects.

    Parameters
    –––-------
        ct_scan: str
            Path to Nifti image file of CT scan.

        brainmask_ct: str
            Path to Nifti image file of corresponding brain mask in CT voxels.

    Returns
    -------
        ct_img: NiBabel image object
            NiBabel image object of CT scan input.

        bm_ct_img: NiBabel image object
            NiBabel image object of brain mask in CT.
    """

    ct_img, bm_ct_img = None, None

    if ct_scan:
        ct_img = nb.load(ct_scan)

    if brainmask_ct:
        bm_ct_img = nb.load(brainmask_ct)

    return ct_img, bm_ct_img


def load_elecs_data(elecfile):
    """
    Load each brain image scan as a NiBabel image object.

    Parameters
    ----------
        elecfile: str
            Path to space-delimited text file of contact labels and contact
            coordinates in mm space.

    Returns
    -------
        eleccoord_mm: dict(str: ndarray)
            Dictionary of contact coordinates in mm space. Maps contact labels
            to contact coordinates, stored as 1x3 numpy arrays.
    """

    eleccoords_mm = {}

    if elecfile.endswith(".txt"):
        with open(elecfile) as f:
            for l in f:
                row = l.split()
                if len(row) == 4:
                    eleccoords_mm[row[0]] = np.array(list(map(float, row[1:])))
                elif len(row) == 6:
                    eleccoords_mm[row[1]] = np.array(list(map(float, row[2:5])))
                else:
                    raise ValueError(
                        "Unrecognized electrode coordinate text format"
                    )
    else:
        matreader = MatReader()
        data = matreader.loadmat(elecfile)

        eleclabels = data["eleclabels"]
        elecmatrix = data["elecmatrix"]
        print(f"Electrode matrix shape: {elecmatrix.shape}")

        for i in range(len(eleclabels)):
            eleccoords_mm[eleclabels[i][0].strip()] = elecmatrix[i]

    print(f'Electrode labels: {eleccoords_mm.keys()}')

    return eleccoords_mm


def apply_atlas(fspatdir, destrieuxfilepath, dktfilepath, fs_lut_fpath):
    """
    Map centroids to an atlas (e.g. Desikan-Killiany, Destriuex) and apply
    white matter and brain masks to label centroids as white matter or out of
    the brain.

    Parameters
    –---------
        fspatdir: str
            Path to freesurfer directory.

        destrieuxfilepath: str
            Path to destrieux atlas for patient.

        dktfilepath: str
            Path to Desikan-Killiany atlas for patient.

        fs_lut_fpath: str
            Path to fs_lut file.

    Returns
    -------
        elec_labels_destriuex: dict(str: ndarray)
            array of contacts labeled with Destriuex atlas.

        elec_labels_DKT: dict(str: ndarray)
            array of contacts labeled with Desikan-Killiany atlas.
    """
    destriuexname = os.path.splitext(os.path.basename(destrieuxfilepath))[0]
    dktname = os.path.splitext(os.path.basename(dktfilepath))[0]

    patid = os.path.basename(os.path.normpath(fspatdir))

    # Apply Atlases, white matter mask, and brainmask
    convert_fsmesh2mlab(subj_dir=os.path.abspath(os.path.dirname(fspatdir)), subj=patid)
    elec_labels_destriuex = label_elecs(subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
                                        subj=patid, hem="lh",
                                        fs_lut_fpath=fs_lut_fpath,
                                        elecfile_prefix=destriuexname, atlas_depth="destriuex"
                                        )
    elec_labels_DKT = label_elecs(subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
                                  subj=patid, hem="lh", fs_lut_fpath=fs_lut_fpath,
                                  elecfile_prefix=dktname, atlas_depth="desikan-killiany"
                                  )
    return elec_labels_destriuex, elec_labels_DKT


def save_organized_elecdict(elecdict, outputfilepath):
    """
    Save centroids as .mat file with attributes eleclabels, which stores
    channel name, electrode type, and depth/grid/strip/other and elecmatrix,
    which stores the centroid coordinates.

    Parameters
    ----------
        elecdict: dict(str: dict(str: np.array))
            Dictionary of channels grouped by electrode.

        outputfilepath: str
            Filepath to save .mat file with annotated electrode labels.
    """
    eleclabels = []
    elecmatrix = []
    for elec in elecdict.keys():
        for chan in elecdict[elec]:
            label = [[chan.strip()], "stereo", "depth"]
            eleclabels.append(label)
            elecmatrix.append(elecdict[elec][chan])
    mat = {"eleclabels": eleclabels, "elecmatrix": elecmatrix}
    scipy.io.savemat(outputfilepath, mat)


def apply_wm_and_brainmask(final_centroids_xyz, atlasfilepath, wmpath, bmpath):
    """
    Apply white matter and brainmask labels to final centroid output and save
    in .mat files.

    Parameters
    ----------
        final_centroids_xyz: dict(str: dict(str: ndarray))
            Dictionary of predicted centroids in xyz (mm) coordinates.

        atlasfilepath: str
            Path to .txt file to save the xyz coordinates of centroids.

        wmpath: str
            Path to white matter mask file.

        bmpath: str
            Path to brain matter mask file.

    Returns
    -------
        anatomy: ndarray
            Anatomy matrix with columns of coordinates, anatomical label, and
            channel label.
    """
    dat = scipy.io.loadmat(atlasfilepath)
    elecmatrix = dat["elecmatrix"]
    anatomy_orig = dat["anatomy"]
    eleclabels = dat["eleclabels"]

    # Load white matter and brain masks
    wm_img = nb.load(wmpath)
    wm_dat = wm_img.get_data()
    bm_img = nb.load(bmpath)
    bm_dat = bm_img.get_data()

    affine = npl.inv(bm_img.affine)

    wm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    bm_label = np.zeros(anatomy_orig.shape[0], dtype=bool)
    # Add two columns in anatomy to store boolean for whether voxel is
    # white matter or brain matter
    anatomy = np.zeros((anatomy_orig.shape[0], anatomy_orig.shape[1] + 2), dtype=object)
    for i, label in enumerate(anatomy_orig):
        chan = str(label[0][0]).strip()
        for elec in final_centroids_xyz:
            if chan in final_centroids_xyz[elec].keys():
                pt = apply_affine(affine, final_centroids_xyz[elec][chan])
                wm_label[i] = wm_dat[list(map(int, pt))] > 0
                bm_label[i] = bm_dat[list(map(int, pt))] > 0
    anatomy[:, :anatomy_orig.shape[1]] = anatomy_orig
    anatomy[:, anatomy_orig.shape[1]] = wm_label
    anatomy[:, anatomy_orig.shape[1] + 1] = bm_label

    save_dict = {"elecmatrix": elecmatrix,
                 "anatomy": anatomy,
                 "eleclabels": eleclabels
                 }
    scipy.io.savemat(atlasfilepath, mdict=save_dict)
    return anatomy


def main(ctimgfile, brainmaskfile, elecinitfile):
    # Hyperparameters
    radius = 4  # radius (in CT voxels) of cylindrical boundary
    threshold = 0.630  # Between 0 and 1. Zeroes voxels with value < threshold
    gap_tolerance = 13  # maximum distance between two adjacent nodes

    # Load data
    print(f"LOADING ELECS DATA FROM: {elecinitfile}")
    elec_coords_mm = load_elecs_data(elecinitfile)
    ct_img, bm_ct_img = load_data(ctimgfile, brainmaskfile)

    # Define pipeline objects to run algorithm
    maskpipe = MaskVolume()
    clusterpipe = Cluster()
    grouppipe = CylindricalGroup()
    postprocesspipe = PostProcessor()

    # Apply masking
    if bm_ct_img:
        brainmasked_ct_img = maskpipe.apply_mask(ct_img, bm_ct_img)
        brainmasked_ct_data = brainmasked_ct_img.get_fdata()
        # Filtering out electrodes not within brainmask
        elecvoxels_in_brain = maskpipe.mask_electrodes(
            elec_coords_mm, brainmasked_ct_img
        )
    else:
        brainmasked_ct_data = ct_img.get_fdata()
        elecvoxels_in_brain = maskpipe.mask_electrodes(
            elec_coords_mm, ct_img
        )
    ct_affine = ct_img.affine

    # Get all voxel clouds per electrode in sorted order
    grouped_contacts_vox = maskpipe.group_contacts(elecvoxels_in_brain)

    # Runs threshold-based clustering algorithm over brainmasked CT
    # for thresholds between 0.63 and 0.64 with a step of 0.005
    clusters, numobj = np.array(
        clusterpipe.find_clusters(brainmasked_ct_data, threshold=threshold)
    )

    # Cluster by cylindrical boundaries
    bound_clusters, sparse_labeled_contacts = grouppipe.cylinder_filter(
        grouped_contacts_vox, clusters, radius
    )

    # Begin postprocessing steps
    processed_clusters = postprocesspipe.process_abnormal_clusters(
        bound_clusters,
        sparse_labeled_contacts
    )

    # Compute centroids for filling gaps
    final_clusters = postprocesspipe.fill_gaps(
        processed_clusters,
        gap_tolerance,
        sparse_labeled_contacts
    )

    # Assign channel labels to the clusters found
    final_labeled_clusters = postprocesspipe.assign_labels(
        final_clusters,
        sparse_labeled_contacts
    )

    # Correct egregious errors with linearly interpolated points
    final_labeled_clusters = postprocesspipe.bruteforce_correction(
        final_labeled_clusters,
        sparse_labeled_contacts
    )

    # Compute centroids for each cluster
    final_centroids_vox = postprocesspipe.cluster_2_centroids(
        final_labeled_clusters
    )

    # Convert final voxels to xyz coordinates
    final_centroids_xyz = postprocesspipe.vox_2_xyz(
        final_centroids_vox,
        ct_affine
    )

    return (
        final_centroids_vox,
        final_centroids_xyz,
        brainmasked_ct_img,
        elecvoxels_in_brain,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ct_nifti_img", help="The CT image volume in its original space."
    )
    parser.add_argument(
        "brainmask_native_file", help="Brain mask mapped to the CT image space."
    )
    parser.add_argument(
        "electrode_initialization_file",
        help="The electrode file with contacts localized to 2 points.",
    )
    # parser.add_argument(
    #     "chanxyz_file",
    #     help="The output datafile with all the electrode centroid points labeled.",
    # )
    parser.add_argument(
        "clustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument("clustered_voxels_file", help="the voxels output datafile")
    parser.add_argument(
        "orgclustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument("orgclustered_voxels_file", help="the voxels output datafile")

    parser.add_argument("binarized_ct_volume", help="The binarized CT volume.")
    parser.add_argument("fsdir", help="The freesurfer output diretroy.")
    parser.add_argument("patid")
    parser.add_argument("fs_lut_fpath")
    parser.add_argument("--wm_native_file", default=None)
    args = parser.parse_args()

    # Extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    brainmask_native_file = args.brainmask_native_file
    electrode_initialization_file = args.electrode_initialization_file
    # chanxyz_file = args.chanxyz_file
    clustered_points_file = args.clustered_points_file
    clustered_voxels_file = args.clustered_voxels_file
    orgclustered_points_file = args.orgclustered_points_file
    orgclustered_voxels_file = args.orgclustered_voxels_file
    binarized_ct_file = args.binarized_ct_volume
    fsdir = args.fsdir
    patid = args.patid
    fs_lut_fpath = args.fs_lut_fpath
    wm_native_file = args.wm_native_file

    # Create electrodes directory if not exist
    elecs_dir = os.path.join(fsdir, "elecs")
    if not os.path.exists(elecs_dir):
        os.mkdir(elecs_dir)

    # Compute the final centroid voxels, centroid xyzs and the binarized CT image volume
    final_centroids_voxels, final_centroids_xyz, binarized_ct_img, _ = main(
        ct_nifti_img, brainmask_native_file, electrode_initialization_file
    )

    # Save output files
    print(f"Saving clustered xyz coords to: {clustered_points_file}.")
    print(f"Saving clustered voxels to: {clustered_voxels_file}.")
    # Save centroids as .mat file with attributes eleclabels, which stores
    save_organized_elecdict(final_centroids_xyz, clustered_points_file)
    save_organized_elecdict(final_centroids_voxels, clustered_voxels_file)
    print(f"Saving binarized CT image volume to: {binarized_ct_file}.")
    binarized_ct_img.to_filename(binarized_ct_file)

    # Save output clustered points to destrieux and dkt atlas files to be labeled
    destrieuxfilepath = os.path.join(
        elecs_dir, "%s_clustered_elec_xyz_destriuex.mat" % (patid)
    )
    dktfilepath = os.path.join(elecs_dir, "%s_clustered_elec_xyz_DK.mat" % (patid))
    save_organized_elecdict(final_centroids_xyz, destrieuxfilepath)
    save_organized_elecdict(final_centroids_xyz, dktfilepath)

    # Output labeled .mat files with atlas, white matter, and brainmask information
    elec_labels_destriuex, elec_labels_DKT = apply_atlas(
        fsdir, destrieuxfilepath, dktfilepath, fs_lut_fpath
    )

    """ SAVE CLUSTERED VOXELS AND POINTS AS TXT FILES WITH CHANNELS PER ROW """
    scipy.io.savemat(orgclustered_voxels_file, mdict={"data": final_centroids_voxels})
    scipy.io.savemat(orgclustered_points_file, mdict={"data": final_centroids_xyz})
