import argparse
import subprocess
import numpy as np


def transform(coords, src_img, dest_img, transform_mat):
    coords_str = " ".join([str(x) for x in coords])

    # print(coords_str)
    cp = subprocess.run(
        "echo %s | img2imgcoord -mm -src %s -dest %s -xfm %s"
        % (coords_str, src_img, dest_img, transform_mat),
        shell=True,
        stdout=subprocess.PIPE,
    )

    transformed_coords = cp.stdout.decode("ascii").strip().split("\n")[-1]
    # print(transformed_coords)
    return np.array([float(x) for x in transformed_coords.split(" ") if x])


def read_label_coords(elecfile):
    labels = []
    labelsxyz = []

    with open(elecfile, "r") as f:
        for _row in f:
            row = _row.split(" ")
            labels.append(row[0])
            labelsxyz.append([float(x) for x in row[1:]])

    return labels, labelsxyz


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("conformed_mri_nifti_img", help="The conformed MRI nifti img.")
    parser.add_argument("mri_nifti_img", help="Brain MRI image space.")
    parser.add_argument(
        "register_conformed_file", help="The conformed mapping transformation file."
    )
    parser.add_argument(
        "clustered_points_file",
        help="The output datafile with all the electrode points clustered.",
    )
    parser.add_argument(
        "outputcoordsfile",
        help="The output datafile for electrodes mapped to correct coords.",
    )
    args = parser.parse_args()

    # extract arguments from parser
    conformed_mri_nifti_img = args.conformed_mri_nifti_img
    mri_nifti_img = args.mri_nifti_img
    register_conformed_file = args.register_conformed_file
    clustered_points_file = args.clustered_points_file
    outputcoordsfile = args.outputcoordsfile

    # read in electrodes file
    labels, labelsxyz = read_label_coords(clustered_points_file)

    # run coordinat transformation
    modified_coords = np.array(
        [
            transform(
                coords, mri_nifti_img, conformed_mri_nifti_img, register_conformed_file
            )
            for coords in labelsxyz
        ]
    )

    with open(outputcoordsfile, "w") as f:
        for i, name in enumerate(labels):
            f.write(
                "%s %.6f %.6f %.6f\n"
                % (
                    name,
                    modified_coords[i][0],
                    modified_coords[i][1],
                    modified_coords[i][2],
                )
            )
