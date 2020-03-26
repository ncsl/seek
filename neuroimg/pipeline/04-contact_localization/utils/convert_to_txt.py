import argparse

import scipy.io


def loadmat(filename):
    """Load .mat file into RAM."""

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
    """Convert coords to 4 x C array."""
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
    clustered_points_file = args.clustered_points_file
    outputcoordsfile = args.outputcoordsfile

    # read in electrodes file
    electxt = read_label_coords(clustered_points_file)

    # write the output to a txt file
    with open(outputcoordsfile, "w") as f:
        for i, name in enumerate(electxt.keys()):
            f.write(
                "%s %.6f %.6f %.6f\n"
                % (name, electxt[name][0], electxt[name][1], electxt[name][2])
            )
