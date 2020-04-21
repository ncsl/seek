import argparse
from pathlib import Path
import scipy.io

from mne_bids.tsv_handler import _from_tsv


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
    parser.add_argument('channels_tsv_fpath')
    parser.add_argument('elecs_dir')
    parser.add_argument(
        "outputcoordsfile",
        help="The output datafile for electrodes mapped to correct coords.",
    )
    args = parser.parse_args()

    # extract arguments from parser
    channels_tsv_fpath = args.channels_tsv_fpath
    elecs_dir = args.elecs_dir
    # clustered_points_file = args.clustered_points_file
    outputcoordsfile = args.outputcoordsfile

    matfiles = [x for x in Path(elecs_dir).glob("*.mat")]
    print(elecs_dir)
    print(matfiles)
    if len(matfiles) > 1:
        raise RuntimeError("There should only be one .mat file inside the "
                           "FreeSurfer elecs directory. This should be created "
                           "from manual annotation of at least entry/exit points "
                           "on each electrode. (Possibly using Fieldtrip Toolbox)")
    clustered_points_file = matfiles[0]

    # read in electrodes file
    electxt = read_label_coords(clustered_points_file)

    print("Saving electrode coordinates as a txt file!")
    print(electxt)

    _electxt = dict()
    channels_tsv = _from_tsv(channels_tsv_fpath)
    for ch in channels_tsv['name']:
        if ch.upper() in electxt:
            _electxt[ch] = electxt[ch]
        else:
            _electxt[ch] = ["n/a", "n/a", "n/a"]
    electxt = _electxt

    # write the output to a txt file
    with open(outputcoordsfile, "w") as f:
        for i, name in enumerate(electxt.keys()):
            f.write(
                # "%s %.6f %.6f %.6f\n"
                "%s %s %s %s\n"
                % (name, electxt[name][0], electxt[name][1], electxt[name][2])
            )
