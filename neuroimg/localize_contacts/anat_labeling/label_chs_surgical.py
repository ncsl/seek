from typing import List, Union

import numpy as np
import nibabel as nb
import scipy.io
import os

from mne_bids import make_bids_basename, make_bids_folders
from mne_bids.tsv_handler import _from_tsv
from nibabel.affines import apply_affine
from nibabel.orientations import aff2axcodes

class MatReader:
    """
    Object to read mat files into a nested dictionary if need be.
    Helps keep strucutre from matlab similar to what is used in python.
    """

    def __init__(self, filename=None):
        self.filename = filename

    def loadmat(self, filename):
        """
        this function should be called instead of direct spio.loadmat
        as it cures the problem of not properly recovering python dictionaries
        from mat files. It calls the function check keys to cure all entries
        which are still mat-objects
        """
        # data = scipy.io.loadmat(
        #     filename,
        #     struct_as_record=False,
        #     squeeze_me=True)

        # load in setup.mat via scipy
        try:
            data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
        except:
            #     setup_dict = mat73.loadmat(_get_setup_fname(source_fpath))
            # finally:
            data = read_matlab(filename)

        return self._check_keys(data)

    def _check_keys(self, dict):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in dict:
            if isinstance(dict[key], scipy.io.matlab.mio5_params.mat_struct):
                dict[key] = self._todict(dict[key])
        return dict

    def _todict(self, matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        dict = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                dict[strg] = self._todict(elem)
            elif isinstance(elem, np.ndarray):
                dict[strg] = self._tolist(elem)
            else:
                dict[strg] = elem
        return dict

    def _tolist(self, ndarray):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, scipy.io.matlab.mio5_params.mat_struct):
                elem_list.append(self._todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(self._tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

def _apply_segmentation_mask(img, seg_mask_arr: np.ndarray) -> nb.Nifti2Image:
    img_data = img.get_fdata()

    # mask image array
    segmented_img_data = np.multiply(img_data, seg_mask_arr)

    # create Nifti image
    segmented_img = nb.Nifti2Image(segmented_img_data, img.affine, header=img.header)
    return segmented_img


def _get_surgical_contacts(
    img: nb.Nifti2Image, ch_names: List, ch_coords: Union[List, np.ndarray]
) -> List:
    """Map mm xyz coordinates into voxel space."""
    affine = img.affine
    inv_affine = np.linalg.inv(affine)
    img_data = img.get_fdata()

    # store surgical channels as a list
    surgical_chs = []

    # map coordinates to voxels and determine if they lay within segmentation volume
    for name, coord in zip(ch_names, ch_coords):
        new_coord = list(map(int, apply_affine(inv_affine, coord)))
        if img_data[new_coord] > 0:
            surgical_chs.append(name)
    return surgical_chs


def _compare_surgical_contacts(surg_chs, clinically_removed_contacts):
    not_found_chs = set(clinically_removed_contacts).difference(set(surg_chs))
    print("Difference: ", not_found_chs)
    return not_found_chs

def read_surgical_mask(fpath):
    fpath = str(fpath)
    if fpath.endswith(".nrrd"):
        surgmask_arr, surgmask_hdr = nrrd.read(fpath, index_order='F')
        # surgmask_arr = np.swapaxes(surgmask_arr, 1, 2)
        space_directions = surgmask_hdr['space directions']
        space_origin = surgmask_hdr['space origin']
        affine = np.vstack((np.hstack((space_directions, space_origin[:, np.newaxis])),
                                np.array([0, 0, 0, 1])))
        # lps_to_ras = np.diag([-1, - 1, 1, 1])
        # ijk_to_ras = lps_to_ras @ ijk_to_lps
        # print("ijk to LPS ", ijk_to_lps)

        # apply to surg mask array
        surgmask_img = nb.Nifti1Image(surgmask_arr, affine)
        print(surgmask_img.affine)
        # get the original orientation
        orig_axcodes = aff2axcodes(surgmask_img.affine)
        orig_ornt = nb.orientations.axcodes2ornt(orig_axcodes,)
                                                 # (('L','R'),('P','A'),('I','S')))

        # get the new orientation and apply it
        # new_ornt = nb.orientations.axcodes2ornt(('L', 'I', 'A'),
        #                                         # (('L','R'),('I','S'),('P','A'))
        #                                         )
        # orig_to_new_ornt = nb.orientations.ornt_transform(orig_ornt, new_ornt)
        #
        # surgmask_arr = nb.apply_orientation(surgmask_arr, orig_to_new_ornt)
        # surgmask_affine = nb.apply_orientation(surgmask_img.affine, orig_to_new_ornt)
        # print("Surgical mask affine: ", surgmask_affine)
        # raise Exception("")
    elif fpath.endswith(".mat"):
        # data_dict = scipy.io.loadmat(fpath)['scirunnrrd']
        loader = MatReader()
        data_dict = loader.loadmat(fpath)['scirunnrrd']
        surgmask_arr = np.array(data_dict['data'])
        surgmask_hdr = data_dict['axis']

    return surgmask_arr, surgmask_hdr

def read_label_coords(elecfile):
    labels = []
    labelsxyz = []

    print("Reading ", elecfile)

    with open(elecfile, "r") as f:
        for _row in f:
            row = _row.split(" ")
            labels.append(row[0])
            labelsxyz.append([float(x) for x in row[1:]])

    return labels, labelsxyz

def read_electrodes_tsv(fname, bids_root, coord_type):
    # read in electrodes
    subj_dir = make_bids_folders(bids_root=bids_root.as_posix(), subject=subject,
                                 session='veeg', make_dir=False)
    fpath = os.path.join(subj_dir, fname)
    print("Reading electrode coord data from: ", fpath)
    electrodes_tsv = _from_tsv(fpath)

    # read in the T1 MRI mapped electrode coordinates
    ch_names = electrodes_tsv['name']
    ch_coords = np.vstack((electrodes_tsv['x'], electrodes_tsv['y'], electrodes_tsv['z'])).T.astype(float)

    # read in the electrode file
    # elecfile = os.path.join(deriv_dir, subject, 'elecs', 'la02_elecxyz_inmri2.txt')
    # ch_names, ch_coords = read_label_coords(elecfile)

    # get channels
    if coord_type == 'vox':
        # convert MRI 'mm' coordinates to voxels
        # ct_ch_vox = apply_affine(np.linalg.inv(ct_affine), ch_coords)
        # ch_voxs = np.array(apply_affine(ct_to_t1_affine, ct_ch_vox)).astype(int)
        ch_voxs = np.array(apply_affine(np.linalg.inv(t1w_affine), ch_coords)).astype(int)
        ch_coords = ch_voxs

    return ch_names, ch_coords

if __name__ == '__main__':
    from pathlib import Path
    import nrrd
    from pprint import pprint
    bids_root = Path("/Users/adam2392/Dropbox/epilepsy_bids/")
    deriv_dir = Path(bids_root / 'derivatives'/ 'freesurfer')
    subject = 'nl04'
    postsurg_dir = Path(deriv_dir / subject / 'postsurgerymri')

    # get the segmented fpath
    segmented_fpath = Path(postsurg_dir / f"{subject}-surgical-segmentation2.nrrd")

    # read the actual mask
    surgmask_arr, surgmask_hdr = read_surgical_mask(segmented_fpath)
    mask_measurement_frame = surgmask_hdr['measurement frame']
    pprint(surgmask_hdr)

    # get the original T1 img
    T1w_fpath = Path(postsurg_dir / "preT1.nii").as_posix()
    t1w_img = nb.load(T1w_fpath)
    t1w_affine = t1w_img.affine
    mm_to_voxel = np.linalg.inv(t1w_affine)

    print(aff2axcodes(mask_measurement_frame))
    print(t1w_img.header)
    print(t1w_img.affine)
    # canonical_img = nb.as_closest_canonical(t1w_img)
    # print("Canonical affine: ", canonical_img.affine)
    # print(nb.aff2axcodes(canonical_img.affine))

    ct_fpath = Path(deriv_dir / subject / 'CT' / 'CT.nii').as_posix()
    ct_img = nb.load(ct_fpath)
    ct_affine = ct_img.affine
    ct_to_t1_affine = np.loadtxt(Path(deriv_dir / subject / 'CT' / 'fsl_ct-to-t1_omat.txt'))

    # get the post -> pre T1 affine
    post_to_pre_affine = np.loadtxt(Path(postsurg_dir / 'fsl_postt1-to-t1_omat.txt'))
    pre_to_post_affine = np.linalg.inv(post_to_pre_affine)

    # get electrodes coordinates in voxel format
    electrodes_fname = make_bids_basename(subject=subject,
                                          session='veeg',
                                          processing='manual',
                                          acquisition='seeg',
                                          suffix="electrodes.tsv")
    ch_names, ch_voxs = read_electrodes_tsv(electrodes_fname, bids_root, coord_type='vox')

    print(surgmask_arr.shape)
    print(len(ch_names), ch_voxs.shape)
    surgmask_arr = np.swapaxes(surgmask_arr, 1, 2)
    # surgmask_arr = np.flip(surgmask_arr, axis=0)
    surgmask_arr = np.flip(surgmask_arr, axis=1)
    surgmask_arr = np.flip(surgmask_arr, axis=2)
    ablation_inds = np.argwhere(surgmask_arr > 0)
    # print(ablation_inds)
    print(ablation_inds.shape)
    surgical_chs = []
    for i, (ch_name, ch_vox) in enumerate(zip(ch_names, ch_voxs)):
        # ch_vox[0] = 256 - ch_vox[0]
        # ch_vox[1] = 256 - ch_vox[1]
        # ch_vox[2] = 256 - ch_vox[2]
        # print(ch_vox)

        # print(surgmask_arr[ch_vox[1], ch_vox[0], ch_vox[2]])
        if surgmask_arr[ch_vox[0], ch_vox[1], ch_vox[2]] == 1:
            surgical_chs.append(ch_name)
        if ch_name.upper() in [
            "L'1", "L'2", "L'3", "L'4",
            "R'1", "R'2", "R'3", "R'4", "R'5", "R'6",
            "H'5", "H'6",  "H'7", "H'8", "H'9", "S'6", "S'7", "S'8", "S'9"
        ]:
            # print(ch_vox)
            dists = []
            for inds in ablation_inds:
                # print(inds)
                dists.append(np.linalg.norm(ch_vox - inds))
            print(ch_name, min(dists))

    print(surgical_chs)