import argparse
import sys
import nibabel as nb
import numpy as np
import imp

sys.path.append('../../image_lddmm_tensorflow/')
sys.path.append('../../ndreg')

from image_lddmm_tensorflow import lddmm


def read_img(fpath):
    img = nb.load(fpath)

    # get info about image space
    if '.img' == fpath[-4:]:
        nxI = img.header['dim'][1:4]
        dxI = img.header['pixdim'][1:4]
    elif fpath.endswith('.mgz'):
        nxI = img.header['dims'][0:3]
        dxI = img.header['delta']
    #     print("Image: ", img[0])
    else:
        # I'm only working with analyze for now
        raise ValueError('Only Analyze images supported for now')

    xI = [np.arange(nxi) * dxi - np.mean(np.arange(nxi) * dxi) for nxi, dxi in zip(nxI, dxI)]
    I = img.get_data()

    return I, xI


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mri_nifti_img', help="Brain MRI image space.")
    parser.add_argument('ct_nifti_img', help="The CT image file.")
    parser.add_argument('output_img_filepath', help="The output datafile with all the electrode points clustered.")
    args = parser.parse_args()

    # extract arguments from parser
    mri_nifti_img = args.mri_nifti_img
    ct_nifti_img = args.ct_nifti_img
    output_img_filepath = args.output_img_filepath

    # read image and dimensions
    I, xI = read_img(mri_nifti_img)
    J, xJ = read_img(ct_nifti_img)

    imp.reload(lddmm)  # for debugging only
    sigmaM = np.std(J)
    out = lddmm.lddmm(I, J,  # atlas and target images
                      xI=xI, xJ=xJ,
                      niter=100,  # number of iterations of gradient descent
                      eV=1e-3 * sigmaM ** 2,  # step size for deformation field update
                      sigmaM=sigmaM,  # noise in image (matching weight 1/2/sigmaM**2)
                      sigmaR=1e1,  # noise in deformation (regularization weight 1/2/sigmaR**2)
                      p=2,  # power of smoothing operator, 2 is typical
                      a=2.0,  # length scale of smoothing operator (mm)
                      nMstep=5,
                      sigmaA=sigmaM,  # std for "artifact" (aka missing tissue)
                      CA0=0.0  # artifact constant value
                      )

    # save the resulting output image filepath
    nb.save(out, output_img_filepath)