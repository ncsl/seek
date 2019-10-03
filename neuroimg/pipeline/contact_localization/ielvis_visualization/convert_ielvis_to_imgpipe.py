import argparse
import sys

import nibabel as nb
import scipy.io
from nibabel.affines import apply_affine

sys.path.append("../../")
from neuroimg.base.utils.data_structures_utils import MatReader


def convert_crs_to_ras(eleccoords, t1_img_file):
    t1_img = nb.load(t1_img_file)
    eleccoords_ras = apply_affine(t1_img.affine, eleccoords)
    return eleccoords_ras


def read_label_coords(elecfile):
    matreader = MatReader()
    elecmat = matreader.loadmat(elecfile)
    return elecmat


def save_as_matfile(elecmat, eleccoords_ras, outputcoordsfile):
    elecmatrix_orig = elecmat["elecmatrix"]
    elec_labels_orig = elecmat["anatomy"]
    elecmontage = elecmat["eleclabels"]

    scipy.io.savemat(
        outputcoordsfile,
        {
            "elecmatrix": elecmatrix_orig,
            "anatomy": elec_labels_orig,
            "eleclabels": elecmontage,
            "elecmatrix_ras": eleccoords_ras,
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mri_nifti_img", help="Brain MRI image space.")
    parser.add_argument(
        "elec_coords_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument(
        "outputcoordsfile",
        help="The output datafile for electrodes mapped to correct coords.",
    )
    args = parser.parse_args()

    # extract arguments from parser
    mri_nifti_img = args.mri_nifti_img
    elec_coords_file = args.elec_coords_file
    outputcoordsfile = args.outputcoordsfile

    # read in electrodes file
    elecmat = read_label_coords(elec_coords_file)
    eleccoords = elecmat["elecmatrix"]

    # apply affine transformation to get the RAS coordinates of file
    eleccoords_ras = convert_crs_to_ras(eleccoords, mri_nifti_img)

    # resave the elec-coords w/ an additional struct for the RAS coordinates
    save_as_matfile(elecmat, eleccoords_ras, outputcoordsfile)
