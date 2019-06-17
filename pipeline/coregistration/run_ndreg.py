import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import nibabel as nb

sys.path.append('../../ndreg/')

import ndreg
from ndreg.ndreg import preprocessor, util


def main(mri_nifti_img, ct_nifti_img, outdir):
    print("Trying to load data.")
    # gets the image spacing in mm for mri
    print(mri_nifti_img, ct_nifti_img)
    t1img = nb.load(mri_nifti_img)
    print("Loaded T1!")
    ctimg = nb.load(ct_nifti_img)
    print("Loaded CT!")

    # extract metadata parameters of the images using nibabel
    atlas_spacing = tuple(t1img.header['pixdim'][0:3].astype('float'))
    image_spacing = tuple(ctimg.header['pixdim'][0:3].astype('float'))
    image_orientation = "".join(nb.aff2axcodes(ctimg.affine)) # return 3-str code
    image_modality = 'lavision'

    params = {
        # input image path
        'image_path': ct_nifti_img,
        # voxel spacing is in mm and corresponds to (x, y, z) spacing
        'image_spacing': image_spacing,
        'image_orientation': image_orientation,
        # the modality can be 'lavision' or 'colm'
        'image_modality': image_modality,
        'atlas_spacing': atlas_spacing,
        'atlas_path': mri_nifti_img,
    }

    # use ndreg.util to read in images
    img = util.imgRead(params['image_path'])
    img.SetSpacing(params['image_spacing'])
    atlas = util.imgRead(params['atlas_path'])
    atlas.SetSpacing(params['atlas_spacing'])

    # print("Running preprocessor")
    # print(len(params['atlas_spacing']), len(params['image_spacing']))
    # print(img.GetDimension(), atlas.GetDimension())
    # print(len(atlas.GetSpacing()))
    # runs preprocessor to
    img_p = preprocessor.preprocess_brain(img,
                                          params['atlas_spacing'],
                                          params['image_modality'],
                                          params['image_orientation'])

    print("Running registration")
    '''
    CURRENTLY METAMORPHISIS FILE CAN'T BE FOUND?
    '''
    # print(img_p.GetDimension())
    # runs registration
    atlas_registered = ndreg.register_brain(atlas, img_p, params['atlas_spacing'], outdir)

    return atlas_registered


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mri_nifti_img', help="Brain MRI image space.")
    parser.add_argument('ct_nifti_img', help="The CT image volume in its original space.")
    parser.add_argument('mapping_transformation_file', help="The mapping transformation file.")
    parser.add_argument('ct_to_mri_nifti_img', help="The output datafile for electrodes mapped to correct coords.")
    parser.add_argument('outdir', help="Output data directory to save results.")
    args = parser.parse_args()

    print("Inside ndreg")

    # extract arguments from parser
    ct_nifti_img = args.ct_nifti_img
    mri_nifti_img = args.mri_nifti_img
    mapping_transformation_file = args.mapping_transformation_file
    ct_to_mri_nifti_img = args.ct_to_mri_nifti_img
    outdir = args.outdir

    # perform registration
    atlas_registered_ct, aff_trans_data = main(mri_nifti_img, ct_nifti_img, outdir)

    # save the image volume
    nb.save(atlas_registered_ct, ct_to_mri_nifti_img)

    # save the mapping transformation file
    affine_file = os.path.join(outdir, 'atlas_to_observed_affine.txt')
    with open(affine_file, 'r') as f:
        affine_data = f.read()

        with open(mapping_transformation_file, 'w') as outf:
            outf.write(affine_data)
