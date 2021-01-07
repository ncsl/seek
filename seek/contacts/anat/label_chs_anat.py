import glob
import os

import nibabel as nb
import numpy as np
import scipy.io

from seek.utils import BidsRoot


def _get_indices_to_exclude(label, indices):
    # These are the indices that won't be used for labeling
    # dont_label = ['EOG','ECG','ROC','LOC','EEG','EKG','NaN','EMG','scalpEEG']
    indices.extend(
        [
            i
            for i, x in enumerate(label)
            if (
                "EOG" in x
                or "ECG" in x
                or "ROC" in x
                or "LOC" in x
                or "EEG" in x
                or "EKG" in x
                or "NaN" in x
                or "EMG" in x
                or x == np.nan
                or "scalpEEG" in x
            )
        ]
    )
    return indices


def _label_grid_and_strips(
    elecmatrix, elec_labels, vert_label, cortex_verts, isnotdepth=None
):
    if isnotdepth is None:
        isnotdepth = []

    # Only use electrodes that are grid or strips
    if len(isnotdepth) > 0:
        elecmatrix_new = elecmatrix[isnotdepth, :]
    else:
        elecmatrix_new = elecmatrix

    print("Finding nearest mesh vertex for each electrode")
    vert_inds, nearest_verts = nearest_electrode_vert(cortex_verts, elecmatrix_new)

    ## Now make a dictionary of the label for each electrode
    elec_labels_notdepth = []
    for v in range(len(vert_inds)):
        if vert_inds[v] in vert_label:
            elec_labels_notdepth.append(vert_label[vert_inds[v]].strip())
        else:
            elec_labels_notdepth.append("Unknown")

    if isnotdepth:
        elec_labels[isnotdepth, 3] = elec_labels_notdepth
        elec_labels[
            np.invert(isnotdepth), 3
        ] = ""  # Set these to an empty string instead of None type
    else:
        elec_labels = np.array(elec_labels_notdepth, dtype=np.object)

    return isnotdepth, elec_labels


def _label_depth(elecmatrix, aparc_dat, LUT, coordinate_type="mm"):
    if coordinate_type == "vox":
        # Define the affine transform to go from surface coordinates to volume coordinates (as CRS, which is
        # the slice *number* as x,y,z in the 3D volume. That is, if there are 256 x 256 x 256 voxels, the
        # CRS coordinate will go from 0 to 255.)
        affine = np.array(
            [
                [-1.0, 0.0, 0.0, 128.0],
                [0.0, 0.0, 1.0, -128.0],
                [0.0, -1.0, 0.0, 128.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
    else:
        affine = np.eye(4)

    # create 4D electrode coordinate array to apply affine transform
    elecs_depths = elecmatrix
    intercept = np.ones(len(elecs_depths))
    elecs_ones = np.column_stack((elecs_depths, intercept))

    # Find voxel CRS
    VoxCRS = (
        np.dot(np.linalg.inv(affine), elecs_ones.transpose()).transpose().astype(int)
    )

    # Make dictionary of labels
    LUT = [row.split() for row in LUT]
    lab = {}
    for row in LUT:
        if (
            len(row) > 1 and row[0][0] != "#" and row[0][0] != "\\"
        ):  # Get rid of the comments
            lname = row[1]
            lab[np.int(row[0])] = lname

    # Label the electrodes according to the aseg volume
    nchans = VoxCRS.shape[0]
    anatomy_labels = np.empty((nchans,), dtype=np.object)
    print("Labeling electrodes...")
    print("Aparc data array has shape: ", aparc_dat.shape)
    print(
        "VoxCRS limits in 3D: ",
        np.max(VoxCRS[:, 0]),
        np.max(VoxCRS[:, 1]),
        np.max(VoxCRS[:, 2]),
    )

    # label each channel
    for elec in np.arange(nchans):
        anatomy_labels[elec] = lab[
            aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]]
        ]
        print(
            "E%d, Vox CRS: [%d, %d, %d], Label #%d = %s"
            % (
                elec,
                VoxCRS[elec, 0],
                VoxCRS[elec, 1],
                VoxCRS[elec, 2],
                aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]],
                anatomy_labels[elec],
            )
        )

    return anatomy_labels


def _split_surf_depth_electrodes(elecmatrix, elecmontage):
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
    indices = []
    indices = _get_indices_to_exclude(long_label, indices)
    indices = _get_indices_to_exclude(short_label, indices)
    indices = _get_indices_to_exclude(grid_or_depth, indices)
    indices.extend(np.where(np.isnan(elecmatrix))[0])
    indices = list(set(indices))
    indices_to_use = list(set(range(len(long_label))) - set(indices))

    # Initialize the cell array that we'll store electrode labels in later
    elec_labels_orig = np.empty((len(long_label), 4), dtype=np.object)
    elec_labels_orig[:, 0] = short_label
    elec_labels_orig[:, 1] = long_label
    elec_labels_orig[:, 2] = grid_or_depth
    elec_labels = np.empty((len(indices_to_use), 4), dtype=np.object)
    elecmatrix = elecmatrix[indices_to_use, :]

    # compute labels for short/long and grid/depth
    short_label = [i for j, i in enumerate(short_label) if j not in indices]
    long_label = [i for j, i in enumerate(long_label) if j not in indices]
    grid_or_depth = [i for j, i in enumerate(grid_or_depth) if j not in indices]

    # make eleclabels a Cx3 array
    elec_labels[:, 0] = short_label
    elec_labels[:, 1] = long_label
    elec_labels[:, 2] = grid_or_depth

    # Find the non depth electrodes
    isnotdepth = np.array([r != "depth" for r in grid_or_depth])
    return isnotdepth, elec_labels, elecmatrix


def nearest_electrode_vert(cortex_verts, elecmatrix):
    """
    Find the vertex on a mesh that is closest to the given electrode coordinates.

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
    """
    nchans = elecmatrix.shape[0]
    d = np.zeros((nchans, cortex_verts.shape[0]))

    # Find the distance between each electrode and all possible vertices
    # on the surface mesh
    for chan in np.arange(nchans):
        d[chan, :] = np.sqrt(
            (elecmatrix[chan, 0] - cortex_verts[:, 0]) ** 2
            + (elecmatrix[chan, 1] - cortex_verts[:, 1]) ** 2
            + (elecmatrix[chan, 2] - cortex_verts[:, 2]) ** 2
        )

    # Find the index of the vertex nearest to each electrode
    vert_inds = np.argmin(d, axis=1)
    nearest_verts = cortex_verts[vert_inds, :]

    return vert_inds, nearest_verts


def convert_fsmesh2mlab(subj_dir, subj, mesh_name="pial"):
    """
    Create surface mesh triangle and vertex .mat files.

    If no argument for mesh_name is given, lh.pial and rh.pial
    are converted into lh_pial_trivert.mat and rh_pial_trivert.mat
    in the Meshes directory (for use in python) and *_lh_pial.mat
    and *_rh_pial.mat for use in MATLAB.

    Parameters
    ----------
    mesh_name : {'pial', 'white', 'inflated'}

    """
    surf_dir = os.path.join(subj_dir, subj, "surf")
    mesh_dir = os.path.join(subj_dir, subj, "Meshes")
    hems = ["lh", "rh"]

    if not os.path.isdir(mesh_dir):
        print("Making Meshes Directory")
        # Make the Meshes directory in subj_dir if it does not yet exist
        os.mkdir(mesh_dir)

    # Loop through hemispheres for this mesh, create one .mat file for each
    for h in hems:
        print("Making %s mesh" % (h))
        mesh_surf = os.path.join(surf_dir, h + "." + mesh_name)
        vert, tri = nb.freesurfer.read_geometry(mesh_surf)
        out_file = os.path.join(mesh_dir, "%s_%s_trivert.mat" % (h, mesh_name))
        out_file_struct = os.path.join(mesh_dir, "%s_%s_%s.mat" % (subj, h, mesh_name))
        scipy.io.savemat(out_file, {"tri": tri, "vert": vert})

        cortex = {"tri": tri + 1, "vert": vert}
        scipy.io.savemat(out_file_struct, {"cortex": cortex})

    if mesh_name == "pial":
        pial_surf_file = dict()
        pial_surf_file["lh"] = os.path.join(mesh_dir, "lh_pial_trivert.mat")
        pial_surf_file["rh"] = os.path.join(mesh_dir, "rh_pial_trivert.mat")
        return pial_surf_file
    else:
        return out_file


def label_elecs_ecog(
    electrodes_names,
    electrodes_array,
    subj_dir,
    subj,
    hem,
    fs_lut_fpath,
    elecfile_prefix="TDT_elecs_all",
    atlas_surf="desikan-killiany",
    atlas_depth="destrieux",
    elecs_all=True,
):
    """
    Automatically labels electrodes based on the freesurfer annotation file.

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
    """
    elecs_dir = os.path.join(subj_dir, subj, "elecs")
    mesh_dir = os.path.join(subj_dir, subj, "Meshes")

    if atlas_surf == "desikan-killiany":
        surf_atlas_flag = ""
    elif atlas_surf == "destrieux":
        surf_atlas_flag = "--a2009s"
    else:
        surf_atlas_flag = ""

    print(subj_dir)
    print(
        "Creating labels from the freesurfer annotation file for use in automated electrode labeling"
    )
    gyri_labels_dir = os.path.join(subj_dir, subj, "label", "gyri")
    if not os.path.isdir(gyri_labels_dir):
        os.mkdir(gyri_labels_dir)

    # This version of mri_annotation2label uses the coarse labels from the Desikan-Killiany Atlas, unless
    # atlas_surf is 'destrieux', in which case the more detailed labels are used
    os.system(
        "mri_annotation2label --subject %s --hemi %s --surface pial %s --outdir %s"
        % (subj, hem, surf_atlas_flag, gyri_labels_dir)
    )
    print("Loading electrode matrix")
    elecfile = os.path.join(elecs_dir, elecfile_prefix + ".mat")
    elecmatrix = scipy.io.loadmat(elecfile)["elecmatrix"]

    # Initialize empty variable for indices of grid and strip electrodes
    isnotdepth = []

    # Choose only the surface or grid electrodes (if not using hd_grid.mat)
    if (
        elecfile_prefix == "TDT_elecs_all"
        or elecfile_prefix == "clinical_elecs_all"
        or elecs_all
    ):
        elecmontage = scipy.io.loadmat(elecfile)["eleclabels"]
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
        indices = [
            i
            for i, x in enumerate(long_label)
            if (
                "EOG" in x
                or "ECG" in x
                or "ROC" in x
                or "LOC" in x
                or "EEG" in x
                or "EKG" in x
                or "NaN" in x
                or "EMG" in x
                or x == np.nan
                or "scalpEEG" in x
            )
        ]
        indices.extend(
            [
                i
                for i, x in enumerate(short_label)
                if (
                    "EOG" in x
                    or "ECG" in x
                    or "ROC" in x
                    or "LOC" in x
                    or "EEG" in x
                    or "EKG" in x
                    or "NaN" in x
                    or "EMG" in x
                    or x == np.nan
                    or "scalpEEG" in x
                )
            ]
        )
        indices.extend(
            [
                i
                for i, x in enumerate(grid_or_depth)
                if (
                    "EOG" in x
                    or "ECG" in x
                    or "ROC" in x
                    or "LOC" in x
                    or "EEG" in x
                    or "EKG" in x
                    or "NaN" in x
                    or "EMG" in x
                    or x == np.nan
                    or "scalpEEG" in x
                )
            ]
        )
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

        short_label_orig, long_label_orig, grid_or_depth_orig = (
            short_label,
            long_label,
            grid_or_depth,
        )
        short_label = [i for j, i in enumerate(short_label) if j not in indices]
        long_label = [i for j, i in enumerate(long_label) if j not in indices]
        grid_or_depth = [i for j, i in enumerate(grid_or_depth) if j not in indices]
        elec_labels[:, 0] = short_label
        elec_labels[:, 1] = long_label
        elec_labels[:, 2] = grid_or_depth

        # Find the non depth electrodes
        isnotdepth = np.array([r != "depth" for r in grid_or_depth])

    # Use the surface label files to get which label goes with each surface vertex
    label_files = glob.glob(os.path.join(gyri_labels_dir, "%s.*.label" % (hem)))
    vert_label = {}
    for label in label_files:
        label_name = label.split(".")[1]
        print("Loading label %s" % label_name)
        fid = open(label, "r")
        d = np.genfromtxt(fid, delimiter=" ", skip_header=2)
        vertnum, x, y, z, junk = d[~np.isnan(d)].reshape((-1, 5)).T
        for v in vertnum:
            vert_label[np.int(v)] = label_name.strip()
        fid.close()

    trivert_file = os.path.join(mesh_dir, "%s_pial_trivert.mat" % (hem))
    cortex_verts = scipy.io.loadmat(trivert_file)["vert"]

    # Only use electrodes that are grid or strips
    if len(isnotdepth) > 0:
        elecmatrix_new = elecmatrix[isnotdepth, :]
    else:
        elecmatrix_new = elecmatrix

    print("Finding nearest mesh vertex for each electrode")
    vert_inds, nearest_verts = nearest_electrode_vert(cortex_verts, elecmatrix_new)

    ## Now make a dictionary of the label for each electrode
    elec_labels_notdepth = []
    for v in range(len(vert_inds)):
        if vert_inds[v] in vert_label:
            elec_labels_notdepth.append(vert_label[vert_inds[v]].strip())
        else:
            elec_labels_notdepth.append("Unknown")

    if (
        elecfile_prefix == "TDT_elecs_all"
        or elecfile_prefix == "clinical_elecs_all"
        or elecs_all
    ):
        elec_labels[isnotdepth, 3] = elec_labels_notdepth
        elec_labels[
            np.invert(isnotdepth), 3
        ] = ""  # Set these to an empty string instead of None type
    else:
        elec_labels = np.array(elec_labels_notdepth, dtype=np.object)
    print("Saving electrode labels for surface electrodes to %s" % (elecfile_prefix))
    ## added by BKD so that elec_mat_grid='hd_grid' works. It does not contain elecmontage
    save_dict = {"elecmatrix": elecmatrix, "anatomy": elec_labels}
    if "elecmontage" in locals():
        save_dict["eleclabels"] = elecmontage
    else:
        print("electmontage does not exist")
    # scipy.io.savemat('%s/%s/elecs/%s'%(subj_dir, subj, elecfile_prefix), save_dict)

    if np.any(np.invert(isnotdepth)):  # If there are depth electrodes, run this part
        print("*************************************************")
        print("Now doing the depth electrodes")

        # Get the volume corresponding to the labels from the Destrieux atlas, which is more
        # detailed than Desikan-Killiany (https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
        if atlas_depth == "desikan-killiany":
            depth_atlas_nm = ""
        elif atlas_depth == "destrieux":
            depth_atlas_nm = ".a2009s"
        else:
            depth_atlas_nm = ".a2009s"

        aseg_file = os.path.join(
            subj_dir, subj, "mri", "aparc%s+aseg.mgz" % (depth_atlas_nm)
        )
        dat = nb.freesurfer.load(aseg_file)
        aparc_dat = dat.get_data()

        # Define the affine transform to go from surface coordinates to volume coordinates (as CRS, which is
        # the slice *number* as x,y,z in the 3D volume. That is, if there are 256 x 256 x 256 voxels, the
        # CRS coordinate will go from 0 to 255.)
        affine = np.array(
            [
                [-1.0, 0.0, 0.0, 128.0],
                [0.0, 0.0, 1.0, -128.0],
                [0.0, -1.0, 0.0, 128.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

        elecs_depths = elecmatrix[np.invert(isnotdepth), :]
        intercept = np.ones(len(elecs_depths))
        elecs_ones = np.column_stack((elecs_depths, intercept))

        # Find voxel CRS
        VoxCRS = (
            np.dot(np.linalg.inv(affine), elecs_ones.transpose())
            .transpose()
            .astype(int)
        )
        # Make meshgrid the same size as aparc_dat (only for gaussian blob version), ignore
        # xx, yy, zz = np.mgrid[0:aparc_dat.shape[0], 0:aparc_dat.shape[1], 0:aparc_dat.shape[2]]
        # unique_labels = np.unique(aparc_dat)
        # unique_labels = unique_labels[unique_labels>0]

        # Get the names of these labels using Freesurfer's lookup table (LUT)
        print("Loading lookup table for freesurfer labels")
        fid = open(fs_lut_fpath)
        LUT = fid.readlines()
        fid.close()

        # Make dictionary of labels
        LUT = [row.split() for row in LUT]
        lab = {}
        for row in LUT:
            if (
                len(row) > 1 and row[0][0] != "#" and row[0][0] != "\\"
            ):  # Get rid of the comments
                lname = row[1]
                lab[np.int(row[0])] = lname

        # Label the electrodes according to the aseg volume
        nchans = VoxCRS.shape[0]
        anatomy = np.empty((nchans,), dtype=np.object)
        print("Labeling electrodes...")

        for elec in np.arange(nchans):
            anatomy[elec] = lab[
                aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]]
            ]
            print(
                "E%d, Vox CRS: [%d, %d, %d], Label #%d = %s"
                % (
                    elec,
                    VoxCRS[elec, 0],
                    VoxCRS[elec, 1],
                    VoxCRS[elec, 2],
                    aparc_dat[VoxCRS[elec, 0], VoxCRS[elec, 1], VoxCRS[elec, 2]],
                    anatomy[elec],
                )
            )

        elec_labels[np.invert(isnotdepth), 3] = anatomy

        # make some corrections b/c of NaNs in elecmatrix
    elec_labels_orig[:, 3] = ""
    elec_labels_orig[indices_to_use, 3] = elec_labels[:, 3]

    print("Saving electrode labels to %s" % (elecfile_prefix))
    scipy.io.savemat(
        os.path.join(elecs_dir, elecfile_prefix + ".mat"),
        {
            "elecmatrix": elecmatrix_orig,
            "anatomy": elec_labels_orig,
            "eleclabels": elecmontage,
        },
    )

    return elec_labels


def label_elecs(
    bids_root,
    ch_names,
    elecmatrix,
    subj_dir,
    subj,
    hem,
    fs_lut_fpath,
    atlas_depth="destrieux",
):
    """
    Automatically labels electrodes based on the freesurfer annotation file.

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
    """
    BidsPatient = BidsRoot(subject_id=subj, bids_root=bids_root)
    elecs_dir = BidsPatient.elecs_dir
    mesh_dir = BidsPatient.mesh_dir
    gyri_labels_dir = BidsPatient.gyri_dir
    mri_dir = BidsPatient.mri_dir

    # get the surface atlas suffix
    # surf_atlas_suffix = BidsPatient.get_surface_atlas_suffix(atlas_surf)

    # Get the volume corresponding to the labels from the Destrieux atlas, which is more
    # detailed than Desikan-Killiany (https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
    depth_atlas_suffix = BidsPatient.get_depth_atlas_suffix(atlas_depth)

    print("# of electrode channels: ", len(ch_names))
    print("Electrode matrix: ", elecmatrix.shape)
    print(subj_dir)
    print(
        "Assigning labels from the freesurfer annotation file for use in automated electrode labeling. "
    )
    # Use the surface label files to get which label goes with each surface vertex
    label_files = glob.glob(os.path.join(gyri_labels_dir, "%s.*.label" % (hem)))
    vert_label = {}
    for label in label_files:
        label_name = label.split(".")[1]
        print("Loading label %s" % label_name)
        fid = open(label, "r")
        d = np.genfromtxt(fid, delimiter=" ", skip_header=2)
        vertnum, x, y, z, junk = d[~np.isnan(d)].reshape((-1, 5)).T
        for v in vertnum:
            vert_label[np.int(v)] = label_name.strip()
        fid.close()

    # load in hemispheric pial triangular-vertices mat files
    trivert_file = os.path.join(mesh_dir, "%s_pial_trivert.mat" % (hem))
    if not os.path.exists(trivert_file):
        raise FileNotFoundError(
            "Trivert file .mat file was not created yet. "
            "Please call `convert_fsmesh2mlab` function first."
        )
    cortex_verts = scipy.io.loadmat(trivert_file)["vert"]

    print("Loading electrode matrix")

    # Initialize empty variable for indices of grid and strip electrodes
    isnotdepth = np.empty((0,), dtype=bool)

    # Choose only the surface or grid electrodes (if not using hd_grid.mat)
    # label seeg
    elecmontage = zip(ch_names, "stereo", "depth")
    isnotdepth, elec_labels, elec_matrix_split = _split_surf_depth_electrodes(
        elecmatrix, elecmontage
    )
    # # label grid and strip electrodes
    # elec_labels = _label_grid_and_strips(elecmatrix, names, vert_label, cortex_verts)
    # print("Saving electrode labels for surface electrodes.")
    # print(isnotdepth)
    # print(elec_labels)

    # if np.any(np.invert(isnotdepth)):  # If there are depth electrodes, run this part
    print("*************************************************")
    print("Now doing the depth electrodes")

    # load in ASEG image file
    aseg_file = os.path.join(mri_dir, "aparc%s+aseg.mgz" % (depth_atlas_suffix))
    depth_atlas_img = nb.freesurfer.load(aseg_file)
    aparc_dat = depth_atlas_img.get_data()

    # Get the names of these labels using Freesurfer's lookup table (LUT)
    print("Loading lookup table for freesurfer labels")
    fid = open(fs_lut_fpath)
    LUT = fid.readlines()
    fid.close()

    # TODO: future proof for when we also label grids in same pipeline
    elecmatrix_depth = elecmatrix
    eleclabels_depth = ch_names
    elec_depth_labels = _label_depth(elecmatrix_depth, aparc_dat, LUT)

    assert len(eleclabels_depth) == len(eleclabels_depth)

    return elec_depth_labels
