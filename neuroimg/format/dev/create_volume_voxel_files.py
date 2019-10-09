import os
import sys
import argparse
import nibabel as nb
import numpy as np

from neuroimg.base.objects import GetSurface, RegionIndexMapping, StructuralDataset

parser = argparse.ArgumentParser()
parser.add_argument("dk_file", help="The filepath to the destrieux parcellation")
parser.add_argument("destrieux_file", help="The filepath to the destrieux parcellation")
parser.add_argument("lut_file", help="The lut table file target.")
parser.add_argument("out_dk_file", help="The output dk filepath.")
parser.add_argument("out_destrieux_file", help="The output destrieux filepath.")


def read_file(filepath):
    img = nb.load(filepath)

    return img


def create_voxel_volume_txt(img, outfilepath):
    # create voxel volume python datastruct
    # img_voxel_volume = img

    # filename = os.path.join(label_direc, hemisphere.value + ".aparc.annot")
    # region_mapping, _, _ = nibabel.freesurfer.io.read_annot(filename)

    np.savetxt(outfilepath, img_voxel_volume, fmt="%d")


if __name__ == "__main__":
    dk_file = parser.dk_file
    destrieux_file = parser.destrieux_file
    lut_file = parser.lut_file
    dk_out_file = parser.out_dk_file
    destrieux_out_file = parser.out_destrieux_file

    # read in the images
    dk_img = read_file(dk_file)
    destrieux_img = read_file(destrieux_file)

    # create voxel volume text files for each atlas
    create_voxel_volume_txt(dk_img, dk_out_file)
    create_voxel_volume_txt(destrieux_img, destrieux_out_file)
