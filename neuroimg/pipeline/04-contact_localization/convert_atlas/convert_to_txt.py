import argparse

import scipy.io


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


def loadmat(filename):
    def _check_keys(dict):
        for key in dict:
            if isinstance(dict[key], scipy.io.matlab.mio5_params.mat_struct):
                dict[key] = _todict(dict[key])
        return dict

    def _todict(matobj):
        dict = {}

        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                dict[strg] = _todict(elem)
            else:
                dict[strg] = elem
        return dict

    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def read_label_coords(elecfilemat):
    print("Reading ", elecfilemat)

    elecmat = loadmat(elecfilemat)
    elecxyz = elecmat["elecf"]

    electxt = {
        elecxyz["label"][i]: list(elecxyz["elecpos"][i])
        for i in range(len(elecxyz["label"]))
    }

    return electxt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("talairach_elec_file", help="Brain MRI image space.")
    parser.add_argument(
        "txt_talairach_file", help="The CT image volume in its original space."
    )
    args = parser.parse_args()

    # extract arguments from parser
    txt_talairach_file = args.txt_talairach_file
    labeled_txt_talairach_file = args.labeled_txt_talairach_file

    # read in electrodes file
    labels, labelsxyz = read_label_coords(txt_talairach_file)
