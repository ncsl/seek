import argparse
import glob
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


def nearest_electrode_vert(cortex_verts, elecmatrix):
    ''' Find the vertex on a mesh that is closest to the given electrode
    coordinates.

    Parameters
    ----------
    cortex_verts : array-like
        [nvertices x 3] matrix of vertices on the cortical surface mesh
    elecmatrix : array-like
        [nchans x 3] matrix of 3D electrode coordinates

    Returns
    -------
    vert_inds : array-like
        Array of vertex indices that are closest to each of the
        electrode
    nearest_verts : array-like
        Coordinates for the nearest cortical vertices
    '''

    nchans = elecmatrix.shape[0]
    d = np.zeros((nchans, cortex_verts.shape[0]))

    # Find the distance between each electrode and all possible vertices
    # on the surface mesh
    for chan in np.arange(nchans):
        d[chan, :] = np.sqrt((elecmatrix[chan, 0] - cortex_verts[:, 0]) ** 2 + \
                             (elecmatrix[chan, 1] - cortex_verts[:, 1]) ** 2 + \
                             (elecmatrix[chan, 2] - cortex_verts[:, 2]) ** 2)

    # Find the index of the vertex nearest to each electrode
    vert_inds = np.argmin(d, axis=1)
    nearest_verts = cortex_verts[vert_inds, :]

    return vert_inds, nearest_verts


def convert_fsmesh2mlab(subj_dir, subj, mesh_name='pial'):
    '''Creates surface mesh triangle and vertex .mat files
    If no argument for mesh_name is given, lh.pial and rh.pial
    are converted into lh_pial_trivert.mat and rh_pial_trivert.mat
    in the Meshes directory (for use in python) and *_lh_pial.mat
    and *_rh_pial.mat for use in MATLAB.

    Parameters
    ----------
    mesh_name : {'pial', 'white', 'inflated'}

    '''
    surf_dir = os.path.join(subj_dir, subj, 'surf')
    mesh_dir = os.path.join(subj_dir, subj, 'Meshes')
    hems = ['lh', 'rh']

    if not os.path.isdir(mesh_dir):
        print('Making Meshes Directory')
        # Make the Meshes directory in subj_dir if it does not yet exist
        os.mkdir(mesh_dir)

    # Loop through hemispheres for this mesh, create one .mat file for each
    for h in hems:
        print("Making %s mesh" % (h))
        mesh_surf = os.path.join(surf_dir, h + '.' + mesh_name)
        vert, tri = nb.freesurfer.read_geometry(mesh_surf)
        out_file = os.path.join(mesh_dir, '%s_%s_trivert.mat' % (h, mesh_name))
        out_file_struct = os.path.join(mesh_dir, '%s_%s_%s.mat' % (subj, h, mesh_name))
        scipy.io.savemat(out_file, {'tri': tri, 'vert': vert})

        cortex = {'tri': tri + 1, 'vert': vert}
        scipy.io.savemat(out_file_struct, {'cortex': cortex})

    if mesh_name == 'pial':
        pial_surf_file = dict()
        pial_surf_file['lh'] = os.path.join(mesh_dir, 'lh_pial_trivert.mat')
        pial_surf_file['rh'] = os.path.join(mesh_dir, 'rh_pial_trivert.mat')
        return pial_surf_file
    else:
        return out_file


def label_elecs(subj_dir, subj, hem, fs_dir, elecfile_prefix='TDT_elecs_all', atlas_surf='desikan-killiany',
                atlas_depth='destrieux',
                elecs_all=True):
    ''' Automatically labels electrodes based on the freesurfer annotation file.
    Assumes TDT_elecs_all.mat or clinical_elecs_all.mat files
    Uses both the Desikan-Killiany Atlas and the Destrieux Atlas, as described
    here: https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation

    Parameters
    ----------
    elecfile_prefix : str, optional
        prefix of the .mat file with the electrode coordinates matrix
    atlas_surf : {'desikan-killiany', 'destrieux'}
        The atlas to use for labeling of surface electrodes.
    atlas_depth : {'destrieux', 'desikan-killiany'}
        The atlas to use for labeling of depth electrodes.
    elecs_all : bool
        Label all electrodes

    Returns
    -------
    elec_labels : array-like
        [nchans x 4] matrix of electrode labels. Columns include short name,
        long name, 'grid'/'depth'/'strip' label, and assigned anatomical label.
    '''
    elecs_dir = os.path.join(subj_dir, subj, "elecs")
    mesh_dir = os.path.join(subj_dir, subj, 'Meshes')

    if atlas_surf == 'desikan-killiany':
        surf_atlas_flag = ''
    elif atlas_surf == 'destrieux':
        surf_atlas_flag = '--a2009s'
    else:
        surf_atlas_flag = ''

    print(subj_dir)
    print('Creating labels from the freesurfer annotation file for use in automated electrode labeling')
    gyri_labels_dir = os.path.join(subj_dir, subj, 'label', 'gyri')
    if not os.path.isdir(gyri_labels_dir):
        os.mkdir(gyri_labels_dir)

    # This version of mri_annotation2label uses the coarse labels from the Desikan-Killiany Atlas, unless
    # atlas_surf is 'destrieux', in which case the more detailed labels are used
    os.system('mri_annotation2label --subject %s --hemi %s --surface pial %s --outdir %s' % (
        subj, hem, surf_atlas_flag, gyri_labels_dir))
    print('Loading electrode matrix')
    elecfile = os.path.join(elecs_dir, elecfile_prefix + '.mat')
    elecmatrix = scipy.io.loadmat(elecfile)['elecmatrix']

    # Initialize empty variable for indices of grid and strip electrodes
    isnotdepth = []

    # Choose only the surface or grid electrodes (if not using hd_grid.mat)
    if elecfile_prefix == 'TDT_elecs_all' or elecfile_prefix == 'clinical_elecs_all' or elecs_all:
        elecmontage = scipy.io.loadmat(elecfile)['eleclabels']
        # Make the cell array into something more usable by python
        short_label = []
        long_label = []
        grid_or_depth = []

        for r in elecmontage:
            short_label.append(r[0][0])  # This is the shortened electrode montage label
            long_label.append(r[1][0])  # This is the long form electrode montage label
            grid_or_depth.append(r[2][0])  # This is the label for grid, depth, or strip

        # These are the indices that won't be used for labeling
        # dont_label = ['EOG','ECG','ROC','LOC','EEG','EKG','NaN','EMG','scalpEEG']
        indices = [i for i, x in enumerate(long_label) if (
                'EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x == np.nan or 'scalpEEG' in x)]
        indices.extend([i for i, x in enumerate(short_label) if (
                'EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x == np.nan or 'scalpEEG' in x)])
        indices.extend([i for i, x in enumerate(grid_or_depth) if (
                'EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x == np.nan or 'scalpEEG' in x)])
        indices.extend(np.where(np.isnan(elecmatrix) == True)[0])
        indices = list(set(indices))
        indices_to_use = list(set(range(len(long_label))) - set(indices))

        # Initialize the cell array that we'll store electrode labels in later
        elec_labels_orig = np.empty((len(long_label), 4), dtype=np.object)
        elec_labels_orig[:, 0] = short_label
        elec_labels_orig[:, 1] = long_label
        elec_labels_orig[:, 2] = grid_or_depth
        elec_labels = np.empty((len(indices_to_use), 4), dtype=np.object)
        elecmatrix_orig = elecmatrix
        elecmatrix = elecmatrix[indices_to_use, :]

        short_label_orig, long_label_orig, grid_or_depth_orig = short_label, long_label, grid_or_depth
        short_label = [i for j, i in enumerate(short_label) if j not in indices]
        long_label = [i for j, i in enumerate(long_label) if j not in indices]
        grid_or_depth = [i for j, i in enumerate(grid_or_depth) if j not in indices]
        elec_labels[:, 0] = short_label
        elec_labels[:, 1] = long_label
        elec_labels[:, 2] = grid_or_depth

        # Find the non depth electrodes
        isnotdepth = np.array([r != 'depth' for r in grid_or_depth])

    # Use the surface label files to get which label goes with each surface vertex
    label_files = glob.glob(os.path.join(gyri_labels_dir, '%s.*.label' % (hem)))
    vert_label = {}
    for label in label_files:
        label_name = label.split('.')[1]
        print('Loading label %s' % label_name)
        fid = open(label, 'r')
        d = np.genfromtxt(fid, delimiter=' ', \
                          skip_header=2)
        vertnum, x, y, z, junk = d[~np.isnan(d)].reshape((-1, 5)).T
        for v in vertnum:
            vert_label[np.int(v)] = label_name.strip()
        fid.close()

    trivert_file = os.path.join(mesh_dir, '%s_pial_trivert.mat' % (hem))
    cortex_verts = scipy.io.loadmat(trivert_file)['vert']

    # Only use electrodes that are grid or strips
    if len(isnotdepth) > 0:
        elecmatrix_new = elecmatrix[isnotdepth, :]
    else:
        elecmatrix_new = elecmatrix

    print('Finding nearest mesh vertex for each electrode')
    vert_inds, nearest_verts = nearest_electrode_vert(cortex_verts, elecmatrix_new)

    ## Now make a dictionary of the label for each electrode
    elec_labels_notdepth = []
    for v in range(len(vert_inds)):
        if vert_inds[v] in vert_label:
            elec_labels_notdepth.append(vert_label[vert_inds[v]].strip())
        else:
            elec_labels_notdepth.append('Unknown')

    if elecfile_prefix == 'TDT_elecs_all' or elecfile_prefix == 'clinical_elecs_all' or elecs_all:
        elec_labels[isnotdepth, 3] = elec_labels_notdepth
        elec_labels[np.invert(isnotdepth), 3] = ''  # Set these to an empty string instead of None type
    else:
        elec_labels = np.array(elec_labels_notdepth, dtype=np.object)
    print('Saving electrode labels for surface electrodes to %s' % (elecfile_prefix))
    ## added by BKD so that elec_mat_grid='hd_grid' works. It does not contain elecmontage
    save_dict = {'elecmatrix': elecmatrix, 'anatomy': elec_labels}
    if 'elecmontage' in locals():
        save_dict['eleclabels'] = elecmontage
    else:
        print('electmontage does not exist')
    # scipy.io.savemat('%s/%s/elecs/%s'%(subj_dir, subj, elecfile_prefix), save_dict)

    if np.any(np.invert(isnotdepth)):  # If there are depth electrodes, run this part
        print('*************************************************')
        print('Now doing the depth electrodes')

        # Get the volume corresponding to the labels from the Destrieux atlas, which is more
        # detailed than Desikan-Killiany (https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
        if atlas_depth == 'desikan-killiany':
            depth_atlas_nm = ''
        elif atlas_depth == 'destrieux':
            depth_atlas_nm = '.a2009s'
        else:
            depth_atlas_nm = '.a2009s'

        aseg_file = os.path.join(subj_dir, subj, 'mri', 'aparc%s+aseg.mgz' % (depth_atlas_nm))
        dat = nb.freesurfer.load(aseg_file)
        aparc_dat = dat.get_data()

        # Define the affine transform to go from surface coordinates to volume coordinates (as CRS, which is
        # the slice *number* as x,y,z in the 3D volume. That is, if there are 256 x 256 x 256 voxels, the
        # CRS coordinate will go from 0 to 255.)
        affine = np.array([[-1., 0., 0., 128.],
                           [0., 0., 1., -128.],
                           [0., -1., 0., 128.],
                           [0., 0., 0., 1.]])

        elecs_depths = elecmatrix[np.invert(isnotdepth), :]
        intercept = np.ones(len(elecs_depths))
        elecs_ones = np.column_stack((elecs_depths, intercept))

        # Find voxel CRS
        VoxCRS = np.dot(np.linalg.inv(affine), elecs_ones.transpose()).transpose().astype(int)
        # Make meshgrid the same size as aparc_dat (only for gaussian blob version), ignore
        # xx, yy, zz = np.mgrid[0:aparc_dat.shape[0], 0:aparc_dat.shape[1], 0:aparc_dat.shape[2]]
        # unique_labels = np.unique(aparc_dat)
        # unique_labels = unique_labels[unique_labels>0]

        # Get the names of these labels using Freesurfer's lookup table (LUT)
        print("Loading lookup table for freesurfer labels")
        fid = open(os.path.join(fs_dir, 'FreeSurferColorLUT.txt'))
        LUT = fid.readlines()
        fid.close()

        # Make dictionary of labels
        LUT = [row.split() for row in LUT]
        lab = {}
        for row in LUT:
            if len(row) > 1 and row[0][0] is not '#' and row[0][0] is not '\\':  # Get rid of the comments
                lname = row[1]
                lab[np.int(row[0])] = lname

        # Label the electrodes according to the aseg volume
        nchans = VoxCRS.shape[0]
        anatomy = np.empty((nchans,), dtype=np.object)
        print("Labeling electrodes...")

        for elec in np.arange(nchans):
            anatomy[elec] = lab[aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]]]
            print("E%d, Vox CRS: [%d, %d, %d], Label #%d = %s" % (
                elec, VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2],
                aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]],
                anatomy[elec]))

        elec_labels[np.invert(isnotdepth), 3] = anatomy

        # make some corrections b/c of NaNs in elecmatrix
    elec_labels_orig[:, 3] = ''
    elec_labels_orig[indices_to_use, 3] = elec_labels[:, 3]

    print('Saving electrode labels to %s' % (elecfile_prefix))
    scipy.io.savemat(os.path.join(elecs_dir, elecfile_prefix + '.mat'), {'elecmatrix': elecmatrix_orig,
                                                                         'anatomy': elec_labels_orig,
                                                                         'eleclabels': elecmontage})

    return elec_labels


# try:
#     sys.path.insert(0, "/Users/ChesterHuynh/img_pipe")
#     from img_pipe import img_pipe
# except ImportError as e:
#     print(e, flush=True)
#     raise Exception("No imgpipe. Run pip install on README for img_pipe.")
#

def load_data(ct_scan, brainmask_ct):
    """
    Load each brain image scan as a NiBabel image object.

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

        brainmask_ct: NiBabel image object
            NiBabel image object of brain mask in CT.
    """

    ct_img = nb.load(ct_scan)
    bm_ct_img = nb.load(brainmask_ct)
    return ct_img, bm_ct_img


def load_elecs_data(elecfile):
    """
    Load each brain image scan as a NiBabel image object.

    Parameters
    ----------
        elecfile: str
            Space-delimited text file of contact labels and contact
            coordinates in mm space.

    Returns
    -------
        elecinitfile: dict()
            A dictionary of contact coordinates in mm space. Keys are
            individual contact labels, and values are the corresponding
            coordinates in mm space.
    """

    elec_coords_mm = {}

    if elecfile.endswith(".txt"):
        with open(elecfile) as f:
            for l in f:
                row = l.split()
                elec_coords_mm[row[0]] = np.array(list(map(float, row)))
    else:
        matreader = MatReader()
        data = matreader.loadmat(elecfile)

        eleclabels = data["eleclabels"]
        elecmatrix = data["elecmatrix"]
        print(f"Electrode matrix shape: {elecmatrix.shape}")

        for i in range(len(eleclabels)):
            elec_coords_mm[eleclabels[i][0].strip()] = elecmatrix[i]

    print(f'Electrode labels: {elec_coords_mm.keys()}')

    return elec_coords_mm


def compute_centroids(chanxyzvoxels):
    """
    Function to return the centroids of each channel label given a list
    of voxels per channel label.

    Parameters
    ----------
        chanxyzvoxels: dict()
            dictionary of electrodes and corresponding dictionary of channels
            to centroid xyz coordinates.

    Returns
    -------
        list of centroids.
    """
    centroids = {}
    for channel, voxels in chanxyzvoxels.items():
        voxels = np.array(voxels)
        centroids[channel] = np.mean(voxels, axis=0)
    return centroids


def apply_atlas(fspatdir, destrieuxfilepath, dktfilepath):
    """
    Map centroids to an atlas (e.g. Desikan-Killiany, Destriuex) and apply
    white matter and brain masks to label centroids as white matter or out of
    the brain.

    Parameters
    –---------
        final_centroids_xyz: dict()
            dictionary of electrodes and corresponding dictionary of channels
            to centroid xyz coordinates.

    Returns
    -------
        elec_labels_destriuex: dict()
            array of contacts labeled with Destriuex atlas.
        elec_labels_DKT: dict()
            array of contacts labeled with Desikan-Killiany atlas.
    """
    destriuexname = os.path.splitext(os.path.basename(destrieuxfilepath))[0]
    dktname = os.path.splitext(os.path.basename(dktfilepath))[0]

    patid = os.path.basename(os.path.normpath(fspatdir))

    # Apply Atlases, white matter mask, and brainmask
    freeCoG = img_pipe.freeCoG(
        subj=patid, subj_dir=os.path.abspath(os.path.dirname(fspatdir)), hem="lh"
    )
    convert_fsmesh2mlab(subj_dir=os.path.abspath(os.path.dirname(fspatdir)), subj=patid)
    elec_labels_destriuex = label_elecs(subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
                                        subj=patid, hem="lh",
                                        fs_dir=os.path.dirname(fspatdir),
                                        elecfile_prefix=destriuexname, atlas_depth="destriuex"
                                        )
    elec_labels_DKT = label_elecs(subj_dir=os.path.abspath(os.path.dirname(fspatdir)),
                                  subj=patid, hem="lh",
                                  fs_dir=os.path.dirname(fspatdir),
                                  elecfile_prefix=dktname, atlas_depth="desikan-killiany"
                                  )
    return elec_labels_destriuex, elec_labels_DKT


def save_organized_elecdict(elecdict, outputfilepath):
    # Save centroids as .mat file with attributes eleclabels, which stores
    # channel name, electrode type, and depth/grid/strip/other and elecmatrix,
    # which stores the centroid coordinates
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
    in .mat files

    Parameters
    ----------
    final_centroids_xyz: dict()
        Dictionary of predicted centroids in xyz (mm) coordinates.

    atlasfilepath: str
        Path to .txt file to save the xyz coordinates of centroids.

    wmpath: str
        Path to white matter mask file.

    bmpath: str
        Path to brain matter mask file.

    Returns
    -------
        Anatomy matrix with columns of coordinates, anatomical label,
        and channel label.
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
    # hyperparameters:
    radius = 4  # radius (in CT voxels) of cylindrical boundary
    threshold = 0.630  # Between 0 and 1. Zero-out voxels with intensity < threshold.
    gap_tolerance = 12.5  # maximum distance between two adjacent nodes

    # Load data
    print("LOADING ELECS DATA FROM: ", elecinitfile)
    elec_coords_mm = load_elecs_data(elecinitfile)
    ct_img, bm_ct_img = load_data(ctimgfile, brainmaskfile)
    ct_data = ct_img.get_fdata()
    bm_ct_data = bm_ct_img.get_fdata()

    # Define pipeline objects to run algorithm
    maskpipe = MaskVolume()
    clusterpipe = Cluster()
    grouppipe = CylindricalGroup()
    postprocesspipe = PostProcessor()

    # Apply masking
    brainmasked_ct_data = maskpipe.apply_mask(bm_ct_data, ct_data)
    brainmasked_ct_img = nb.Nifti1Image(brainmasked_ct_data, bm_ct_img.affine)
    ct_affine = ct_img.affine

    # Filtering out electrodes not within brainmask
    elecvoxels_in_brain = maskpipe.filter_electrodes_bm(
        elec_coords_mm, brainmasked_ct_img
    )
    # Get all voxel clouds per electrode
    voxels_per_electrode = maskpipe.sort_contacts(elecvoxels_in_brain)

    # Runs threshold-based clustering algorithm over brainmasked CT
    # for thresholds between 0.62 and 0.65 with a step of 0.005
    clusters, numobj = np.array(clusterpipe.find_clusters(brainmasked_ct_data))

    # Cluster by cylindrical boundaries
    clusters_by_cylinder, sparse_elec_labels, sparse_elec_coords = grouppipe.cylinder_filter(
        elecvoxels_in_brain, clusters[threshold], radius
    )

    # Begin postprocessing steps
    processed_clusters = postprocesspipe.process_abnormal_clusters(
        clusters[threshold],
        elecvoxels_in_brain,
        clusters_by_cylinder,
        sparse_elec_labels,
    )

    # Compute centroids for filling gaps
    centroids = {}
    for electrode in processed_clusters:
        centroid = compute_centroids(processed_clusters[electrode])
        centroids[electrode] = centroid
    centroids = postprocesspipe.reassign_labels(centroids)

    dists = postprocesspipe.compute_dists(centroids)
    final_centroids_voxels, dists = postprocesspipe.fill_gaps(
        centroids, dists, gap_tolerance
    )

    # convert final voxels to xyz coordinates
    final_centroids_xyz = postprocesspipe.vox_2_xyz(final_centroids_voxels, ct_affine)

    return (
        final_centroids_voxels,
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
    parser.add_argument(
        "chanxyz_file",
        help="The output datafile with all the electrode centroid points labeled.",
    )
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
    parser.add_argument("--wm_native_file", default=None)
    args = parser.parse_args()

    # extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    brainmask_native_file = args.brainmask_native_file
    electrode_initialization_file = args.electrode_initialization_file
    chanxyz_file = args.chanxyz_file
    clustered_points_file = args.clustered_points_file
    clustered_voxels_file = args.clustered_voxels_file
    orgclustered_points_file = args.orgclustered_points_file
    orgclustered_voxels_file = args.orgclustered_voxels_file
    binarized_ct_file = args.binarized_ct_volume
    fsdir = args.fsdir
    patid = args.patid
    wm_native_file = args.wm_native_file

    # create electrodes directory if not exist
    elecs_dir = os.path.join(fsdir, "elecs")
    if not os.path.exists(elecs_dir):
        os.mkdir(elecs_dir)

    # compute the final centroid voxels, centroid xyzs and the binarized CT image volume.
    final_centroids_voxels, final_centroids_xyz, binarized_ct_img, elecvoxels_in_brain = main(
        ct_nifti_img, brainmask_native_file, electrode_initialization_file
    )

    # save output files
    print(f"Saving clustered xyz coords to: {clustered_points_file}.")
    print(f"Saving clustered voxels to: {clustered_voxels_file}.")
    # Save centroids as .mat file with attributes eleclabels, which stores
    save_organized_elecdict(final_centroids_xyz, clustered_points_file)
    save_organized_elecdict(final_centroids_voxels, clustered_voxels_file)
    print(f"Saving binarized CT image volume to: {binarized_ct_file}.")
    binarized_ct_img.to_filename(binarized_ct_file)

    # save output clustered points to destrieux and dkt atlas files to be labeled
    destrieuxfilepath = os.path.join(
        elecs_dir, "%s_clustered_elec_xyz_destriuex.mat" % (patid)
    )
    dktfilepath = os.path.join(elecs_dir, "%s_clustered_elec_xyz_DK.mat" % (patid))
    save_organized_elecdict(final_centroids_xyz, destrieuxfilepath)
    save_organized_elecdict(final_centroids_xyz, dktfilepath)

    # Output labeled .mat files with atlas, white matter, and brainmask information
    elec_labels_destriuex, elec_labels_DKT = apply_atlas(
        fsdir, destrieuxfilepath, dktfilepath
    )

    # LOOKS LIKE THIS IS REPEATING BUT FOR WM AND BRAINMASK?
    # apply_wm_and_brainmask(final_centroids_xyz, destrieuxfilepath, wm_native_file, brainmask_native_file)
    # apply_wm_and_brainmask(final_centroids_xyz, dktfilepath, wm_native_file, brainmask_native_file)

    """ SAVE CLUSTERED VOXELS AND POINTS AS TXT FILES WITH CHANNELS PER ROW """
    scipy.io.savemat(orgclustered_voxels_file, mdict={"data": final_centroids_voxels})
    scipy.io.savemat(orgclustered_points_file, mdict={"data": final_centroids_xyz})
    # with open(orgclustered_voxels_file, 'w') as f:
    #     for elec in final_centroids_voxels:
    #         for chan in final_centroids_voxels[elec]:
    #             vox = final_centroids_voxels[elec][chan]
    #             f.write("%s %.6f %.6f %.6f\n" % (chan, vox[0], vox[1], vox[2]))
    # with open(orgclustered_points_file, 'w') as f:
    #     for elec in final_centroids_xyz:
    #         for chan in final_centroids_xyz[elec]:
    #             pt = final_centroids_xyz[elec][chan]
    #             f.write("%s %.6f %.6f %.6f\n" % (chan, pt[0], pt[1], pt[2]))
